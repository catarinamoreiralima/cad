#define _POSIX_C_SOURCE 200809L

#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <omp.h>

enum Estado {
    NAO_COMBUSTIVEL, INTACTA, EM_CHAMAS, QUEIMADA, CONTENCAO, NUM_ESTADOS
};

typedef struct {
    int L, C, P, T, seed, limiar;
    int vento_linha, vento_coluna, intensidade;
    int F, Z;

    int *cobertura;
    int *umidade;
    int *ativacao;
    int *estado_atual;
    int *proximo_estado;
    int *tempo_atual;
    int *proximo_tempo;

    int pesos[8];
} Simulacao;

typedef struct {
    int passos_executados;
    int total_ignicoes;
    int pico_passo;
    int pico_qtd;
} Resultado;

static const int DL[8] = {-1, -1, -1,  0, 0,  1, 1, 1};
static const int DC[8] = {-1,  0,  1, -1, 1, -1, 0, 1};

static int *alocar_vetor(size_t n) {
    int *vetor = n <= SIZE_MAX / sizeof(*vetor)
        ? calloc(n ? n : 1, sizeof(*vetor))
        : NULL;
    if (vetor == NULL) {
        fprintf(stderr, "Erro: falha na alocacao de memoria.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        exit(EXIT_FAILURE);
    }
    return vetor;
}

// Usa a mesma sequência da referência sequencial, inclusive no macOS.
static int fire_rand_r(unsigned int *seed) {
#if defined(__APPLE__)
    unsigned int next = *seed;
    unsigned int resultado;

    next = next * 1103515245U + 12345U;
    resultado = (next / 65536U) % 2048U;
    next = next * 1103515245U + 12345U;
    resultado = (resultado << 10) ^ ((next / 65536U) % 1024U);
    next = next * 1103515245U + 12345U;
    resultado = (resultado << 10) ^ ((next / 65536U) % 1024U);
    *seed = next;
    return (int)resultado;
#else
    return rand_r(seed);
#endif
}

static int duracao_queima(int cobertura) {
    return cobertura == 2 ? 2 : 4;
}

static void gerar_floresta(Simulacao *sim, int *combustiveis_iniciais) {
    int n = sim->L * sim->C;
    unsigned int seed = (unsigned int)sim->seed;
    sim->cobertura = alocar_vetor(n);
    sim->umidade = alocar_vetor(n);
    sim->ativacao = alocar_vetor(n);
    sim->estado_atual = alocar_vetor(n);
    sim->tempo_atual = alocar_vetor(n);

    for (int i = 0; i < n; i++) {
        int valor = fire_rand_r(&seed) % 100;
        int cobertura;
        if (valor < 10) cobertura = 0;
        else if (valor < 20) cobertura = 1;
        else if (valor < 55) cobertura = 2;
        else cobertura = 3;
        sim->cobertura[i] = cobertura;
        sim->umidade[i] = fire_rand_r(&seed) % 101;
        sim->estado_atual[i] = cobertura >= 2 ? INTACTA : NAO_COMBUSTIVEL;
        sim->ativacao[i] = -1;
        *combustiveis_iniciais += cobertura >= 2;
    }
}

// Retorna o diagnóstico da entrada, ou NULL se a floresta estiver pronta.
static const char *ler_entrada(FILE *entrada, Simulacao *sim,
                               int *combustiveis_iniciais) {
    if (fscanf(entrada, "%d%d%d%d%d%d", &sim->L, &sim->C, &sim->P,
               &sim->T, &sim->seed, &sim->limiar) != 6 ||
        sim->L <= 0 || sim->C <= 0 || sim->P < 0 ||
        sim->T <= 0 || sim->limiar <= 0) {
        return "parametros gerais invalidos";
    }

    int L = sim->L;
    int C = sim->C;
    // Scatterv e Gatherv usam contagens e deslocamentos do tipo int.
    if (L > INT_MAX / C) {
        return "dimensoes excedem o limite de contagem MPI";
    }
    if (fscanf(entrada, "%d%d%d", &sim->vento_linha, &sim->vento_coluna,
               &sim->intensidade) != 3 ||
        sim->vento_linha < -1 || sim->vento_linha > 1 ||
        sim->vento_coluna < -1 || sim->vento_coluna > 1 ||
        (sim->vento_linha == 0 && sim->vento_coluna == 0) ||
        sim->intensidade < 0 || sim->intensidade > 5) {
        return "configuracao do vento invalida";
    }

    if (fscanf(entrada, "%d%d", &sim->F, &sim->Z) != 2 ||
        sim->F < 0 || sim->F > L * C || sim->Z < 0) {
        return "quantidade de focos ou zonas invalida";
    }

    gerar_floresta(sim, combustiveis_iniciais);
    for (int f = 0; f < sim->F; f++) {
        int linha, coluna;
        if (fscanf(entrada, "%d%d", &linha, &coluna) != 2 ||
            linha < 0 || linha >= L || coluna < 0 || coluna >= C) {
            return "coordenadas de foco invalidas";
        }
        int idx = linha * C + coluna;
        if (sim->estado_atual[idx] == NAO_COMBUSTIVEL) {
            return "foco sobre celula nao combustivel";
        }
        if (sim->estado_atual[idx] == EM_CHAMAS) {
            return "foco repetido";
        }
        sim->estado_atual[idx] = EM_CHAMAS;
        sim->tempo_atual[idx] = duracao_queima(sim->cobertura[idx]);
    }

    for (int z = 0; z < sim->Z; z++) {
        int passo_ativacao, linha_inicial, coluna_inicial, linha_final, coluna_final;
        if (fscanf(entrada, "%d%d%d%d%d", &passo_ativacao, &linha_inicial,
                   &coluna_inicial, &linha_final, &coluna_final) != 5 ||
            passo_ativacao < 0 || passo_ativacao >= sim->P ||
            linha_inicial < 0 || linha_inicial > linha_final || linha_final >= L ||
            coluna_inicial < 0 || coluna_inicial > coluna_final || coluna_final >= C) {
            return "limites ou passo da zona de contencao invalidos";
        }
        for (int linha = linha_inicial; linha <= linha_final; linha++) {
            for (int coluna = coluna_inicial; coluna <= coluna_final; coluna++) {
                int idx = linha * C + coluna;
                if (sim->ativacao[idx] == -1 || passo_ativacao < sim->ativacao[idx]) {
                    sim->ativacao[idx] = passo_ativacao;
                }
            }
        }
    }
    return NULL;
}

static void precalcular_pesos(Simulacao *sim) {
    for (int k = 0; k < 8; k++) {
        int basico = abs(DL[k]) + abs(DC[k]) == 1 ? 10 : 7;
        int alinhamento = -DL[k] * sim->vento_linha - DC[k] * sim->vento_coluna;
        int peso = basico + sim->intensidade * alinhamento;
        sim->pesos[k] = peso < 1 ? 1 : peso;
    }
}

static int potencial_ignicao(const Simulacao *sim, int idx) {
    int C = sim->C;
    int coluna = idx % C;
    // O índice local da primeira linha própria inclui a linha fantasma superior.
    ptrdiff_t centro = (ptrdiff_t)idx + C;
    int soma = 0;

    #pragma omp simd reduction(+:soma)
    for (int k = 0; k < 8; k++) {
        int vizinha = coluna + DC[k];
        if (vizinha >= 0 && vizinha < C) {
            ptrdiff_t j = centro + (ptrdiff_t)DL[k] * C + DC[k];
            soma += (sim->estado_atual[j] == EM_CHAMAS) * sim->pesos[k];
        }
    }
    int fator = sim->cobertura[idx] == 2 ? 8 : 12;
    return soma * fator * (100 - sim->umidade[idx]) / 100;
}

static void trocar_fronteiras(int *estado, int C, int n, int acima, int abaixo) {
    // MPI_PROC_NULL deixa a linha fantasma da borda externa zerada.
    MPI_Sendrecv(estado + C, C, MPI_INT, acima, 0,
                 estado + (size_t)n + C, C, MPI_INT, abaixo, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(estado + n, C, MPI_INT, abaixo, 1,
                 estado, C, MPI_INT, acima, 1,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

static void reduzir_estatisticas(int ignicoes, int chamas, int globais[2]) {
    int locais[2] = {ignicoes, chamas};
    MPI_Allreduce(locais, globais, 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}

static void trocar_matrizes(Simulacao *sim) {
    int *tmp_estado = sim->estado_atual;
    sim->estado_atual = sim->proximo_estado;
    sim->proximo_estado = tmp_estado;

    int *tmp_tempo = sim->tempo_atual;
    sim->tempo_atual = sim->proximo_tempo;
    sim->proximo_tempo = tmp_tempo;
}

static Resultado simular(Simulacao *sim, int n, int acima, int abaixo) {
    Resultado resultado = {.pico_passo = -1};
    int C = sim->C;
    int continua = sim->P > 0 && sim->F > 0;
    int ignicoes_no_passo = 0;
    int em_chamas_no_passo = 0;

    if (!continua) {
        return resultado;
    }

    #pragma omp parallel num_threads(sim->T) default(none) \
        shared(sim, C, n, acima, abaixo, resultado, continua, ignicoes_no_passo, em_chamas_no_passo)
    {
        while (continua) {
            #pragma omp for schedule(static)
            for (int idx = 0; idx < n; idx++) {
                size_t idx_estado = (size_t)idx + C;
                if (sim->ativacao[idx] == resultado.passos_executados &&
                    sim->estado_atual[idx_estado] == INTACTA) {
                    sim->estado_atual[idx_estado] = CONTENCAO;
                }
            }

            // A mesma thread que inicializou MPI executa toda a comunicação.
            if (omp_get_thread_num() == 0 && n > 0) {
                trocar_fronteiras(sim->estado_atual, C, n, acima, abaixo);
            }
            #pragma omp barrier

            #pragma omp for schedule(guided) reduction(+:ignicoes_no_passo,em_chamas_no_passo)
            for (int idx = 0; idx < n; idx++) {
                size_t idx_estado = (size_t)idx + C;
                int estado = sim->estado_atual[idx_estado];
                int tempo = 0;
                if (estado == EM_CHAMAS) {
                    tempo = sim->tempo_atual[idx] - 1;
                    if (tempo == 0) estado = QUEIMADA;
                } else if (estado == INTACTA &&
                           potencial_ignicao(sim, idx) >= sim->limiar) {
                    estado = EM_CHAMAS;
                    tempo = duracao_queima(sim->cobertura[idx]);
                    ignicoes_no_passo++;
                }
                sim->proximo_estado[idx_estado] = estado;
                sim->proximo_tempo[idx] = tempo;
                em_chamas_no_passo += estado == EM_CHAMAS;
            }

            if (omp_get_thread_num() == 0) {
                int globais[2];
                reduzir_estatisticas(ignicoes_no_passo, em_chamas_no_passo, globais);

                resultado.total_ignicoes += globais[0];
                if (globais[0] > resultado.pico_qtd) {
                    resultado.pico_qtd = globais[0];
                    resultado.pico_passo = resultado.passos_executados;
                }

                trocar_matrizes(sim);
                resultado.passos_executados++;
                continua = resultado.passos_executados < sim->P && globais[1] > 0;
                ignicoes_no_passo = 0;
                em_chamas_no_passo = 0;
            }
            #pragma omp barrier
        }
    }

    return resultado;
}

static void imprimir_resultado(const Simulacao *global, int n, int combustiveis_iniciais,
                               Resultado r, double tempo) {
    int contagem[NUM_ESTADOS] = {0};
    unsigned long long checksum = 0;
    for (int i = 0; i < n; i++) {
        contagem[global->estado_atual[i]]++;
        checksum = checksum * 31ULL + (unsigned long long)global->estado_atual[i];
        checksum = checksum * 31ULL + (unsigned long long)global->tempo_atual[i];
    }

    double queimado = 0.0;
    double protegido = 0.0;
    if (combustiveis_iniciais > 0) {
        queimado = 100.0 * (contagem[QUEIMADA] + contagem[EM_CHAMAS]) / combustiveis_iniciais;
        protegido = 100.0 * contagem[CONTENCAO] / combustiveis_iniciais;
    }

    printf("passos: %d\n", r.passos_executados);
    printf("nao_combustiveis: %d\n", contagem[NAO_COMBUSTIVEL]);
    printf("intactas: %d\n", contagem[INTACTA]);
    printf("em_chamas: %d\n", contagem[EM_CHAMAS]);
    printf("queimadas: %d\n", contagem[QUEIMADA]);
    printf("contencao: %d\n", contagem[CONTENCAO]);
    printf("total_ignicoes: %d\n", r.total_ignicoes);
    printf("pico_ignicoes: %d %d\n", r.pico_passo, r.pico_qtd);
    printf("percentual_queimado: %.2f\n", queimado);
    printf("percentual_protegido: %.2f\n", protegido);
    printf("checksum: %llu\n", checksum);
    printf("tempo: %.6f\n", tempo);
}

static void liberar_estruturas(Simulacao *sim) {
    free(sim->cobertura);
    free(sim->umidade);
    free(sim->ativacao);
    free(sim->estado_atual);
    free(sim->tempo_atual);
    free(sim->proximo_estado);
    free(sim->proximo_tempo);
}

// Empacota somente os inteiros da configuração, sem enviar ponteiros.
static void difundir_parametros(const Simulacao *global, Simulacao *sim) {
    int parametros[] = {
        global->L, global->C, global->P, global->T, global->seed, global->limiar,
        global->vento_linha, global->vento_coluna, global->intensidade, global->F
    };
    MPI_Bcast(parametros, sizeof(parametros) / sizeof(parametros[0]),
              MPI_INT, 0, MPI_COMM_WORLD);

    sim->L = parametros[0];
    sim->C = parametros[1];
    sim->P = parametros[2];
    sim->T = parametros[3];
    sim->seed = parametros[4];
    sim->limiar = parametros[5];
    sim->vento_linha = parametros[6];
    sim->vento_coluna = parametros[7];
    sim->intensidade = parametros[8];
    sim->F = parametros[9];
}

int main(int argc, char *argv[]) {
    int suporte;
    if (MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &suporte) != MPI_SUCCESS) {
        return EXIT_FAILURE;
    }

    int rank, processos;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processos);
    if (suporte < MPI_THREAD_FUNNELED) {
        fprintf(stderr, "Erro: MPI_THREAD_FUNNELED indisponivel.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        return EXIT_FAILURE;
    }

    Simulacao global = {0};
    int combustiveis_iniciais = 0;
    int entrada_valida = 1;
    if (rank == 0) {
        const char *erro = NULL;
        if (argc != 2) {
            erro = "informe um unico arquivo de entrada";
        } else {
            FILE *entrada = fopen(argv[1], "r");
            if (entrada == NULL) {
                erro = "nao foi possivel abrir o arquivo de entrada";
            } else {
                erro = ler_entrada(entrada, &global, &combustiveis_iniciais);
                fclose(entrada);
            }
        }
        if (erro != NULL) {
            fprintf(stderr, "Erro: %s.\n", erro);
            entrada_valida = 0;
        }
    }
    MPI_Bcast(&entrada_valida, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!entrada_valida) {
        liberar_estruturas(&global);
        MPI_Finalize();
        return EXIT_FAILURE;
    }
    Simulacao sim = {0};
    difundir_parametros(&global, &sim);

    int L = sim.L;
    int C = sim.C;
    int linhas_base = L / processos;
    int resto = L % processos;
    int n = (linhas_base + (rank < resto)) * C;
    int ativos = L < processos ? L : processos;
    int acima = rank > 0 && rank < ativos ? rank - 1 : MPI_PROC_NULL;
    int abaixo = rank < ativos - 1 ? rank + 1 : MPI_PROC_NULL;
    int halo = n > 0 ? C : 0;
    size_t tamanho_estado = (size_t)n + 2 * (size_t)halo;

    int *contagens = NULL;
    int *deslocamentos = NULL;
    if (rank == 0) {
        contagens = alocar_vetor(processos);
        deslocamentos = alocar_vetor(processos);
        int inicio = 0;
        for (int r = 0; r < processos; r++) {
            contagens[r] = (linhas_base + (r < resto)) * C;
            deslocamentos[r] = inicio;
            inicio += contagens[r];
        }
    }

    sim.cobertura = alocar_vetor(n);
    sim.umidade = alocar_vetor(n);
    sim.ativacao = alocar_vetor(n);
    sim.estado_atual = alocar_vetor(tamanho_estado);
    sim.proximo_estado = alocar_vetor(tamanho_estado);
    sim.tempo_atual = alocar_vetor(n);
    sim.proximo_tempo = alocar_vetor(n);
    int *origens[] = {global.cobertura, global.umidade, global.ativacao,
                      global.estado_atual, global.tempo_atual};
    int *destinos[] = {sim.cobertura, sim.umidade, sim.ativacao,
                       sim.estado_atual + halo, sim.tempo_atual};
    for (size_t i = 0; i < sizeof(origens) / sizeof(origens[0]); i++) {
        MPI_Scatterv(origens[i], contagens, deslocamentos, MPI_INT,
                     destinos[i], n, MPI_INT, 0, MPI_COMM_WORLD);
    }
    free(global.cobertura);
    free(global.umidade);
    free(global.ativacao);
    global.cobertura = NULL;
    global.umidade = NULL;
    global.ativacao = NULL;

    precalcular_pesos(&sim);
    omp_set_dynamic(0);

    MPI_Barrier(MPI_COMM_WORLD);
    double inicio = omp_get_wtime();
    Resultado resultado = simular(&sim, n, acima, abaixo);
    double duracao = omp_get_wtime() - inicio;

    double tempo;
    MPI_Reduce(&duracao, &tempo, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Gatherv(sim.estado_atual + halo, n, MPI_INT, global.estado_atual, contagens,
                deslocamentos, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gatherv(sim.tempo_atual, n, MPI_INT, global.tempo_atual, contagens,
                deslocamentos, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        imprimir_resultado(&global, L * C, combustiveis_iniciais, resultado, tempo);
    }

    liberar_estruturas(&sim);
    liberar_estruturas(&global);
    free(contagens);
    free(deslocamentos);
    MPI_Finalize();
    return EXIT_SUCCESS;
}
