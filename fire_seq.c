#if !defined(_WIN32) && !defined(_POSIX_C_SOURCE)
#define _POSIX_C_SOURCE 199506L
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <omp.h> // usado apenas para medir tempo com omp_get_wtime() (secao 12 do enunciado)

// rand_r e' POSIX e nao existe no MinGW/Windows; no macOS existe, mas usa um
// algoritmo diferente do da glibc e gera outra floresta para a mesma seed.
// Nesses dois casos usa-se esta versao de compatibilidade, que reproduz o
// algoritmo do rand_r da glibc, mantendo a mesma floresta e checksum do Linux.
// Nos demais sistemas POSIX, usa-se a funcao nativa.
static int fire_rand_r(unsigned int *seed) {
#if defined(_WIN32) || defined(__APPLE__)
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

// Representa um foco inicial de incendio (celula que comeca em chamas).
typedef struct {
    int linha;
    int coluna;
} Foco;

// Representa uma zona de contencao: um retangulo (limites inclusivos) que,
// no inicio do passo "passo_ativacao", transforma celulas intactas em contencao.
typedef struct {
    int passo_ativacao;
    int linha_inicial;
    int coluna_inicial;
    int linha_final;
    int coluna_final;
} ZonaContencao;

// Estrutura central da simulacao: parametros de entrada, focos/zonas e as
// matrizes (armazenadas linearmente, indice = linha * C + coluna).
typedef struct {
    int L, C, P, T, seed, limiar;
    int vento_linha, vento_coluna, intensidade;

    // quantidade de focos iniciais e zonas de contenção
    int F, Z;
    Foco *focos;
    ZonaContencao *zonas;

    // ponteiros para as matrizes principais da floresta
    int *cobertura;       // tipo de cobertura de cada celula (fixo): 0 agua, 1 solo exposto, 2 rasteira, 3 floresta
    int *umidade;         // umidade de cada celula (fixa), 0 a 100
    int *estado_atual;    // estado da celula no passo atual (ver Quadro 6.4.1: 0 nao combustivel, 1 intacta, 2 em chamas, 3 queimada, 4 contencao)
    int *proximo_estado;  // estado calculado para o proximo passo, a partir do estado_atual
    int *tempo_atual;     // tempo de queima restante da celula no passo atual
    int *proximo_tempo;   // tempo de queima calculado para o proximo passo
    int *ativacao;        // ativacao[idx] = -1 (fora de zona) ou o passo em que a zona ativa essa celula

    int pesos[8];          // peso pv de cada uma das 8 direcoes de Moore, pre-calculado uma unica vez
} Simulacao;

// Deslocamentos fixos dos 8 vizinhos de Moore (ordem usada por potencial_ignicao/precalcular_pesos).
static const int DL[8] = {-1, -1, -1,  0, 0,  1, 1, 1};
static const int DC[8] = {-1,  0,  1, -1, 1, -1, 0, 1};

// Le e valida o arquivo de entrada (secoes 4 e 5 do enunciado).
// Encerra o programa com EXIT_FAILURE em qualquer inconsistencia.
void ler_entrada(const char *arquivo_entrada, Simulacao *sim) {
    FILE *arquivo = fopen(arquivo_entrada, "r");
    if (arquivo == NULL) {
        fprintf(stderr, "Erro: Nao foi possivel abrir o arquivo %s\n", arquivo_entrada);
        exit(EXIT_FAILURE);
    }

    // 1a linha: configuracao geral (secao 4.1)
    if (fscanf(arquivo, "%d %d %d %d %d %d",
               &sim->L, &sim->C, &sim->P, &sim->T, &sim->seed, &sim->limiar) != 6) {
        fprintf(stderr, "Erro na leitura da primeira linha.\n");
        exit(EXIT_FAILURE);
    }

    if (sim->L <= 0 || sim->C <= 0 || sim->P < 0 || sim->T <= 0 || sim->limiar <= 0) {
        fprintf(stderr, "Erro: Parametros gerais invalidos.\n");
        exit(EXIT_FAILURE);
    }

    // Todos os indices da simulacao sao int. Rejeitar uma matriz cujo numero
    // de celulas nao cabe nesse tipo evita overflow em L*C e nos indices.
    long long total_celulas = (long long)sim->L * (long long)sim->C;
    if (total_celulas > INT_MAX || (size_t)total_celulas > SIZE_MAX / sizeof(int)) {
        fprintf(stderr, "Erro: Dimensoes da matriz excedem o limite suportado.\n");
        exit(EXIT_FAILURE);
    }

    // 2a linha: configuracao do vento (secao 4.2)
    if (fscanf(arquivo, "%d %d %d",
               &sim->vento_linha,
               &sim->vento_coluna,
               &sim->intensidade) != 3)
    {
        fprintf(stderr, "Erro: configuracao do vento incompleta ou invalida.\n");
        fclose(arquivo);
        exit(EXIT_FAILURE);
    }
    if (sim->vento_linha < -1 || sim->vento_linha > 1 ||
        sim->vento_coluna < -1 || sim->vento_coluna > 1 ||
        (sim->vento_linha == 0 && sim->vento_coluna == 0)) {
        fprintf(stderr, "Erro: Direcao do vento invalida.\n");
        exit(EXIT_FAILURE);
    }
    if (sim->intensidade < 0 || sim->intensidade > 5) {
        fprintf(stderr, "Erro: Intensidade do vento fora do limite [0-5].\n");
        exit(EXIT_FAILURE);
    }

    // 3a linha: quantidade de focos (F) e zonas (Z) (secao 4.3)
    if (fscanf(arquivo, "%d %d", &sim->F, &sim->Z) != 2)
    {
        fprintf(stderr, "Erro: quantidades de focos e zonas invalidas.\n");
        fclose(arquivo);
        exit(EXIT_FAILURE);
    }
    if (sim->F < 0 || sim->Z < 0) {
        fprintf(stderr, "Erro: Quantidade de focos ou zonas nao pode ser negativa.\n");
        exit(EXIT_FAILURE);
    }

    if ((size_t)sim->F > SIZE_MAX / sizeof(Foco) ||
        (size_t)sim->Z > SIZE_MAX / sizeof(ZonaContencao)) {
        fprintf(stderr, "Erro: Quantidade de focos ou zonas excede o limite suportado.\n");
        exit(EXIT_FAILURE);
    }

    sim->focos = malloc((size_t)sim->F * sizeof(*sim->focos));
    sim->zonas = malloc((size_t)sim->Z * sizeof(*sim->zonas));
    if ((sim->F > 0 && sim->focos == NULL) || (sim->Z > 0 && sim->zonas == NULL)) {
        fprintf(stderr, "Erro: Falha na alocacao de focos ou zonas.\n");
        free(sim->focos);
        free(sim->zonas);
        exit(EXIT_FAILURE);
    }

    // proximas F linhas: focos iniciais de incendio
    for (int i = 0; i < sim->F; i++) {
        if (fscanf(arquivo, "%d %d",
                   &sim->focos[i].linha,
                   &sim->focos[i].coluna) != 2)
        {
            fprintf(stderr, "Erro: foco %d incompleto ou invalido.\n", i + 1);
            fclose(arquivo);
            free(sim->focos);
            free(sim->zonas);
            exit(EXIT_FAILURE);
        }

        if (sim->focos[i].linha < 0 || sim->focos[i].linha >= sim->L ||
            sim->focos[i].coluna < 0 || sim->focos[i].coluna >= sim->C) {
            fprintf(stderr, "Erro: Foco fora da matriz.\n");
            exit(EXIT_FAILURE);
        }

        // verifica se o foco ja foi informado antes (foco repetido)
        for (int j = 0; j < i; j++) {
            if (sim->focos[i].linha == sim->focos[j].linha &&
                sim->focos[i].coluna == sim->focos[j].coluna) {
                fprintf(stderr, "Erro: Foco repetido encontrado.\n");
                exit(EXIT_FAILURE);
            }
        }
    }

    // proximas Z linhas: zonas de contencao
    for (int i = 0; i < sim->Z; i++) {
        if (fscanf(arquivo, "%d %d %d %d %d",
                   &sim->zonas[i].passo_ativacao,
                   &sim->zonas[i].linha_inicial,
                   &sim->zonas[i].coluna_inicial,
                   &sim->zonas[i].linha_final,
                   &sim->zonas[i].coluna_final) != 5)
        {
            fprintf(stderr, "Erro: zona %d incompleta ou invalida.\n", i + 1);
            fclose(arquivo);
            free(sim->focos);
            free(sim->zonas);
            exit(EXIT_FAILURE);
        }

        if (sim->zonas[i].passo_ativacao < 0 || sim->zonas[i].passo_ativacao >= sim->P) {
            fprintf(stderr, "Erro: Passo de ativacao da zona invalido.\n");
            exit(EXIT_FAILURE);
        }

        if (sim->zonas[i].linha_inicial > sim->zonas[i].linha_final ||
            sim->zonas[i].coluna_inicial > sim->zonas[i].coluna_final) {
            fprintf(stderr, "Erro: Limites iniciais da zona superiores aos finais.\n");
            exit(EXIT_FAILURE);
        }

        if (sim->zonas[i].linha_inicial < 0 || sim->zonas[i].linha_final >= sim->L ||
            sim->zonas[i].coluna_inicial < 0 || sim->zonas[i].coluna_final >= sim->C) {
            fprintf(stderr, "Erro: Zona de contencao fora da matriz.\n");
            exit(EXIT_FAILURE);
        }
    }

    fclose(arquivo);
}

// Aloca todos os vetores de tamanho L*C usados pela simulacao.
// Nao entra no trecho cronometrado (secao 12).
void alocar_estruturas(Simulacao *sim) {
    size_t n = (size_t)sim->L * (size_t)sim->C;

    sim->cobertura = malloc(n * sizeof(*sim->cobertura));
    sim->umidade = malloc(n * sizeof(*sim->umidade));
    sim->estado_atual = malloc(n * sizeof(*sim->estado_atual));
    sim->proximo_estado = malloc(n * sizeof(*sim->proximo_estado));
    sim->tempo_atual = malloc(n * sizeof(*sim->tempo_atual));
    sim->proximo_tempo = malloc(n * sizeof(*sim->proximo_tempo));
    sim->ativacao = malloc(n * sizeof(*sim->ativacao));

    if (sim->cobertura == NULL || sim->umidade == NULL ||
        sim->estado_atual == NULL || sim->proximo_estado == NULL ||
        sim->tempo_atual == NULL || sim->proximo_tempo == NULL ||
        sim->ativacao == NULL) {
        fprintf(stderr, "Erro: Falha na alocacao de memoria.\n");
        exit(EXIT_FAILURE);
    }
}

// Gera a cobertura e a umidade de cada celula (secoes 6.2 e 6.3).
// Percorre a matriz em ordem crescente de linha e coluna com uma unica
// sequencia de rand_r, pois o resultado precisa ser deterministico e
// identico entre as versoes sequencial e paralela.
void gerar_floresta(Simulacao *sim) {
    // rand_r exige unsigned int*. Usar uma variavel do tipo correto evita a
    // violacao de aliasing causada por converter int* para unsigned int*.
    unsigned int seed = (unsigned int)sim->seed;

    for (int linha = 0; linha < sim->L; linha++) {
        for (int coluna = 0; coluna < sim->C; coluna++) {
            int idx = linha * sim->C + coluna;

            // cobertura: 0-9 agua, 10-19 solo exposto, 20-54 rasteira, 55-99 floresta (Quadro 6.2.1)
            int valor = fire_rand_r(&seed) % 100;
            if (valor < 10) sim->cobertura[idx] = 0;
            else if (valor < 20) sim->cobertura[idx] = 1;
            else if (valor < 55) sim->cobertura[idx] = 2;
            else sim->cobertura[idx] = 3;

            // umidade gerada logo em seguida, ainda na mesma celula (ordem exigida pelo enunciado)
            sim->umidade[idx] = fire_rand_r(&seed) % 101;

            // agua/solo exposto = nao combustivel (0); rasteira/floresta comecam intactas (1)
            sim->estado_atual[idx] = (sim->cobertura[idx] <= 1) ? 0 : 1;
            sim->tempo_atual[idx] = 0;
        }
    }

}

// Aplica os focos iniciais de incendio (secao 6.5): marca as celulas como
// "em chamas" com o tempo de queima inicial do seu tipo de cobertura.
void aplicar_focos(Simulacao *sim) {
    for (int i = 0; i < sim->F; i++) {
        int idx = sim->focos[i].linha * sim->C + sim->focos[i].coluna;

        // foco sobre agua/solo exposto e' entrada invalida (secao 6.5)
        if (sim->cobertura[idx] == 0 || sim->cobertura[idx] == 1) {
            fprintf(stderr, "Erro: Foco de incendio sobre celula nao combustivel.\n");
            exit(EXIT_FAILURE);
        }

        sim->estado_atual[idx] = 2; // em chamas
        sim->tempo_atual[idx] = (sim->cobertura[idx] == 2) ? 2 : 4; // 2 passos p/ rasteira, 4 p/ floresta
    }
}

// Constroi o vetor ativacao[idx] (secao 6.6): -1 se a celula nao pertence a
// nenhuma zona, ou o passo de ativacao caso contrario. Em zonas sobrepostas,
// prevalece o menor passo de ativacao.
void construir_mapa_ativacao(Simulacao *sim) {
    long long n = (long long)sim->L * (long long)sim->C;
    for (long long i = 0; i < n; i++) {
        sim->ativacao[i] = -1;
    }

    for (int z = 0; z < sim->Z; z++) {
        ZonaContencao *zona = &sim->zonas[z];
        for (int linha = zona->linha_inicial; linha <= zona->linha_final; linha++) {
            for (int coluna = zona->coluna_inicial; coluna <= zona->coluna_final; coluna++) {
                long long idx = (long long)linha * sim->C + coluna;
                if (sim->ativacao[idx] == -1 || zona->passo_ativacao < sim->ativacao[idx]) {
                    sim->ativacao[idx] = zona->passo_ativacao;
                }
            }
        }
    }
}

// Pre-calcula o peso pv de cada uma das 8 direcoes de Moore (secao 8): esses
// valores dependem apenas de vento_linha/vento_coluna/intensidade, fixos
// durante toda a simulacao, entao calcula-los uma unica vez evita repetir a
// mesma conta a cada vizinho, de cada celula intacta, em cada passo.
// Chamada em main() logo apos ler_entrada, fora do trecho cronometrado.
void precalcular_pesos(Simulacao *sim) {
    for (int k = 0; k < 8; k++) {
        int prop_linha = -DL[k];
        int prop_coluna = -DC[k];
        int ortogonal = (abs(prop_linha) + abs(prop_coluna) == 1);
        int p_basico = ortogonal ? 10 : 7;

        int A = prop_linha * sim->vento_linha + prop_coluna * sim->vento_coluna;
        int pv = p_basico + sim->intensidade * A;
        if (pv < 1) pv = 1;

        sim->pesos[k] = pv;
    }
}

// Ativa, no inicio do passo "passo", as zonas de contencao programadas para
// esse passo (secao 7.2 / Quadro 7.2.1). Apenas celulas intactas viram
// contencao; os demais estados (inclusive "em chamas") permanecem como estao,
// pois a contencao nao pode apagar um incendio ja iniciado.
void ativar_zonas(Simulacao *sim, int passo) {
    long long n = (long long)sim->L * (long long)sim->C;
    for (long long i = 0; i < n; i++) {
        if (sim->ativacao[i] == passo && sim->estado_atual[i] == 1) {
            sim->estado_atual[i] = 4;
        }
    }
}

// Calcula o potencial de ignicao I de uma celula intacta (secao 8), somando
// a contribuicao de cada vizinho de Moore em chamas conforme a direcao do
// vento e o tipo de posicao (ortogonal/diagonal). Os pesos pv de cada
// direcao ja vem pre-calculados em sim->pesos (ver precalcular_pesos).
int potencial_ignicao(Simulacao *sim, long long idx) {
    long long linha = idx / sim->C;
    long long coluna = idx % sim->C;
    int S = 0;

    // percorre os 8 vizinhos de Moore, ignorando os que ficam fora da matriz
    for (int k = 0; k < 8; k++) {
        long long vl = linha + DL[k];
        long long vc = coluna + DC[k];
        if (vl < 0 || vl >= sim->L || vc < 0 || vc >= sim->C) continue;

        long long vidx = vl * sim->C + vc;
        if (sim->estado_atual[vidx] != 2) continue; // so vizinhos em chamas contribuem

        S += sim->pesos[k];
    }

    // fator de combustivel: 8 para vegetacao rasteira, 12 para floresta (Quadro 6.2.2)
    int fator = (sim->cobertura[idx] == 2) ? 8 : 12;
    return (S * fator * (100 - sim->umidade[idx])) / 100; // divisao inteira (floor)
}

// Calcula o proximo_estado/proximo_tempo de todas as celulas a partir do
// estado_atual (secao 7.3), contando quantas novas ignicoes ocorreram nesse
// passo. Nunca le nem escreve estado_atual/tempo_atual, apenas os proximos.
void calcular_proximo_estado(Simulacao *sim, int *ignicoes_no_passo) {
    long long n = (long long)sim->L * (long long)sim->C;
    *ignicoes_no_passo = 0;

    for (long long idx = 0; idx < n; idx++) {
        int estado = sim->estado_atual[idx];

        if (estado == 0 || estado == 3 || estado == 4) {
            // nao combustivel, queimada e contencao permanecem como estao
            sim->proximo_estado[idx] = estado;
            sim->proximo_tempo[idx] = 0;
        } else if (estado == 2) {
            // em chamas: reduz o tempo de queima; ao chegar a zero, vira queimada
            int novo_tempo = sim->tempo_atual[idx] - 1;
            if (novo_tempo == 0) {
                sim->proximo_estado[idx] = 3;
                sim->proximo_tempo[idx] = 0;
            } else {
                sim->proximo_estado[idx] = 2;
                sim->proximo_tempo[idx] = novo_tempo;
            }
        } else {
            // intacta: entra em chamas se o potencial de ignicao atingir o limiar
            int I = potencial_ignicao(sim, idx);
            if (I >= sim->limiar) {
                sim->proximo_estado[idx] = 2;
                sim->proximo_tempo[idx] = (sim->cobertura[idx] == 2) ? 2 : 4;
                (*ignicoes_no_passo)++;
            } else {
                sim->proximo_estado[idx] = 1;
                sim->proximo_tempo[idx] = 0;
            }
        }
    }
}

// Condicao de parada (secao 9): verdadeiro enquanto existir alguma celula em chamas.
int existe_em_chamas(Simulacao *sim) {
    long long n = (long long)sim->L * (long long)sim->C;
    for (long long i = 0; i < n; i++) {
        if (sim->estado_atual[i] == 2) return 1;
    }
    return 0;
}

// Troca os papeis de atual/proximo (passo 4 da secao 7.1) apenas trocando os
// ponteiros, sem copiar os vetores.
void trocar_matrizes(Simulacao *sim) {
    int *tmp_estado = sim->estado_atual;
    sim->estado_atual = sim->proximo_estado;
    sim->proximo_estado = tmp_estado;

    int *tmp_tempo = sim->tempo_atual;
    sim->tempo_atual = sim->proximo_tempo;
    sim->proximo_tempo = tmp_tempo;
}

// Conta as celulas inicialmente combustiveis (vegetacao rasteira + floresta),
// usada como base para os percentuais queimado/protegido (secao 10).
int contar_combustiveis_iniciais(Simulacao *sim) {
    long long n = (long long)sim->L * (long long)sim->C;
    int total = 0;
    for (long long i = 0; i < n; i++) {
        if (sim->cobertura[i] == 2 || sim->cobertura[i] == 3) total++;
    }
    return total;
}

// Conta quantas celulas terminaram a simulacao em cada um dos 5 estados.
void contar_estados_finais(Simulacao *sim, int *nao_combustiveis, int *intactas,
                            int *em_chamas, int *queimadas, int *contencao) {
    long long n = (long long)sim->L * (long long)sim->C;
    *nao_combustiveis = 0;
    *intactas = 0;
    *em_chamas = 0;
    *queimadas = 0;
    *contencao = 0;

    for (long long i = 0; i < n; i++) {
        switch (sim->estado_atual[i]) {
            case 0: (*nao_combustiveis)++; break;
            case 1: (*intactas)++; break;
            case 2: (*em_chamas)++; break;
            case 3: (*queimadas)++; break;
            case 4: (*contencao)++; break;
        }
    }
}

// Calcula o checksum final (secao 10), sequencialmente e na ordem linear da
// matriz, exatamente como especificado no enunciado.
unsigned long long calcular_checksum(Simulacao *sim) {
    unsigned long long checksum = 0;
    long long n = (long long)sim->L * (long long)sim->C;

    for (long long i = 0; i < n; i++) {
        checksum = checksum * 31ULL + (unsigned long long)sim->estado_atual[i];
        checksum = checksum * 31ULL + (unsigned long long)sim->tempo_atual[i];
    }

    return checksum;
}

// Imprime o resultado final exatamente no formato exigido pela secao 11.
void imprimir_resultado(int passos, int nao_combustiveis, int intactas, int em_chamas,
                         int queimadas, int contencao, int total_ignicoes,
                         int pico_passo, int pico_qtd, double percentual_queimado,
                         double percentual_protegido, unsigned long long checksum,
                         double tempo) {
    printf("passos: %d\n", passos);
    printf("nao_combustiveis: %d\n", nao_combustiveis);
    printf("intactas: %d\n", intactas);
    printf("em_chamas: %d\n", em_chamas);
    printf("queimadas: %d\n", queimadas);
    printf("contencao: %d\n", contencao);
    printf("total_ignicoes: %d\n", total_ignicoes);
    printf("pico_ignicoes: %d %d\n", pico_passo, pico_qtd);
    printf("percentual_queimado: %.2f\n", percentual_queimado);
    printf("percentual_protegido: %.2f\n", percentual_protegido);
    printf("checksum: %llu\n", checksum);
    printf("tempo: %.6f\n", tempo);
}

// Libera toda a memoria alocada em ler_entrada e alocar_estruturas.
void liberar_estruturas(Simulacao *sim) {
    free(sim->focos);
    free(sim->zonas);
    free(sim->cobertura);
    free(sim->umidade);
    free(sim->estado_atual);
    free(sim->proximo_estado);
    free(sim->tempo_atual);
    free(sim->proximo_tempo);
    free(sim->ativacao);
}

int main(int argc, char *argv[]){

    if (argc != 2) {
        printf("Erro de passagem de argumento.\n Uso: ./fire_seq entrada.txt");
        exit(EXIT_FAILURE);
    }
    Simulacao sim;

    // preparacao da simulacao: nada disto entra no trecho cronometrado (secao 12)
    ler_entrada(argv[1], &sim);
    precalcular_pesos(&sim);
    alocar_estruturas(&sim);
    gerar_floresta(&sim);
    aplicar_focos(&sim);
    construir_mapa_ativacao(&sim);

    int combustiveis_iniciais = contar_combustiveis_iniciais(&sim);

    int passo = 0;
    int total_ignicoes = 0;
    int pico_passo = -1; // -1 indica "nenhuma ignicao ocorreu" (secao 10)
    int pico_qtd = 0;

    double t_inicio = omp_get_wtime();

    // se nao houver nenhuma celula em chamas logo apos a inicializacao,
    // nenhum passo e' executado (secao 9)
    if (existe_em_chamas(&sim)) {
        while (passo < sim.P) {
            // ordem de execucao de cada passo, conforme secao 7.1:
            ativar_zonas(&sim, passo);          // 1. ativar zonas programadas para este passo

            int ignicoes_no_passo;
            calcular_proximo_estado(&sim, &ignicoes_no_passo); // 2. calcular proximo estado

            // 3. estatisticas do proximo estado (total e pico de ignicoes; empate mantem o 1o passo)
            total_ignicoes += ignicoes_no_passo;
            if (ignicoes_no_passo > pico_qtd) {
                pico_qtd = ignicoes_no_passo;
                pico_passo = passo;
            }

            trocar_matrizes(&sim); // 4. trocar as matrizes
            passo++;

            // 5. verificar condicao de parada
            if (!existe_em_chamas(&sim)) {
                break;
            }
        }
    }

    double t_fim = omp_get_wtime();
    int passos_executados = passo;

    // calculo dos resultados finais: fora do trecho cronometrado (secao 12)
    int nao_combustiveis, intactas, em_chamas, queimadas, contencao;
    contar_estados_finais(&sim, &nao_combustiveis, &intactas, &em_chamas, &queimadas, &contencao);

    double percentual_queimado = 0.0;
    double percentual_protegido = 0.0;
    if (combustiveis_iniciais > 0) {
        // celulas em chamas contam como queimadas pois ja foram atingidas pelo incendio (secao 10)
        percentual_queimado = 100.0 * (queimadas + em_chamas) / combustiveis_iniciais;
        percentual_protegido = 100.0 * contencao / combustiveis_iniciais;
    }

    unsigned long long checksum = calcular_checksum(&sim);

    imprimir_resultado(passos_executados, nao_combustiveis, intactas, em_chamas, queimadas,
                        contencao, total_ignicoes, pico_passo, pico_qtd, percentual_queimado,
                        percentual_protegido, checksum, t_fim - t_inicio);

    liberar_estruturas(&sim);
    return 0;
}