#if !defined(_WIN32) && !defined(_POSIX_C_SOURCE)
#define _POSIX_C_SOURCE 199506L
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <omp.h>

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
// Identica a versao sequencial: leitura/validacao nao e' paralelizada.
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
// Permanece sequencial mesmo na versao paralela: o enunciado exige
// explicitamente que essa geracao seja feita por uma unica thread, e ela
// nao faz parte do trecho cronometrado.
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

// Calcula o potencial de ignicao I de uma celula intacta (secao 8).
// Reescrito de forma "branchless" (sem continue) sobre os 8 vizinhos de
// Moore pre-calculados, para permitir vetorizacao com omp simd.
int potencial_ignicao(Simulacao *sim, long long idx) {
    long long linha = idx / sim->C;
    long long coluna = idx % sim->C;
    int S = 0;

    #pragma omp simd reduction(+:S)
    for (int k = 0; k < 8; k++) {
        long long vl = linha + DL[k];
        long long vc = coluna + DC[k];
        int dentro = (vl >= 0 && vl < sim->L && vc >= 0 && vc < sim->C);
        // indice "seguro": se o vizinho estiver fora da matriz, reaproveita
        // idx (sempre valido) em vez de acessar fora dos limites; o fator
        // 'dentro' abaixo garante que a contribuicao desse vizinho seja zero
        long long vidx = dentro ? (vl * sim->C + vc) : idx;
        int em_chamas = dentro && (sim->estado_atual[vidx] == 2);

        S += em_chamas * sim->pesos[k];
    }

    int fator = (sim->cobertura[idx] == 2) ? 8 : 12;
    return (S * fator * (100 - sim->umidade[idx])) / 100;
}

// Calcula o proximo_estado/proximo_tempo de UMA celula a partir do
// estado_atual (secao 7.3), sem nenhuma diretiva omp: a paralelizacao fica
// por conta de quem chama (o laco "omp for" em simular()), ja que cada
// celula so escreve na sua propria posicao (sem dependencia entre indices).
// Retorna 1 se houve uma nova ignicao nesta celula; "ficou_em_chamas" indica
// se a celula estara em chamas no proximo passo.
int atualizar_celula(Simulacao *sim, long long idx, int *ficou_em_chamas) {
    int estado = sim->estado_atual[idx];
    int ignizou = 0;
    *ficou_em_chamas = 0;

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
            *ficou_em_chamas = 1;
        }
    } else {
        // intacta: entra em chamas se o potencial de ignicao atingir o limiar
        int I = potencial_ignicao(sim, idx);
        if (I >= sim->limiar) {
            sim->proximo_estado[idx] = 2;
            sim->proximo_tempo[idx] = (sim->cobertura[idx] == 2) ? 2 : 4;
            ignizou = 1;
            *ficou_em_chamas = 1;
        } else {
            sim->proximo_estado[idx] = 1;
            sim->proximo_tempo[idx] = 0;
        }
    }

    return ignizou;
}

// Condicao de parada (secao 9), usada apenas para a checagem inicial (se ha
// chamas logo apos a inicializacao); dentro do laco, "tem_chamas" ja vem
// calculado por reducao em calcular_proximo_estado.
int existe_em_chamas(Simulacao *sim) {
    long long n = (long long)sim->L * (long long)sim->C;
    for (long long i = 0; i < n; i++) {
        if (sim->estado_atual[i] == 2) return 1;
    }
    return 0;
}

// Troca os papeis de atual/proximo (passo 4 da secao 7.1) apenas trocando os
// ponteiros. Executada por uma unica thread (dentro de um single), entao
// nenhuma thread le sim->estado_atual/tempo_atual enquanto a troca acontece.
void trocar_matrizes(Simulacao *sim) {
    int *tmp_estado = sim->estado_atual;
    sim->estado_atual = sim->proximo_estado;
    sim->proximo_estado = tmp_estado;

    int *tmp_tempo = sim->tempo_atual;
    sim->tempo_atual = sim->proximo_tempo;
    sim->proximo_tempo = tmp_tempo;
}

// Regiao paralela persistente: abre uma unica vez (evita o custo de criar e
// destruir threads a cada passo) e executa todo o laco da simulacao dentro
// dela. A ordem de execucao de cada passo segue a secao 7.1: 1. ativar
// zonas; 2. calcular proximo estado; 3. estatisticas; 4. trocar matrizes;
// 5. verificar condicao de parada. A ativacao de zona (passo 1) e' inlineada
// no topo do corpo do "omp for" principal (passo 2): como ela so transiciona
// INTACTA(1)->CONTENCAO(4) e nunca toca em EM_CHAMAS(2), a ordem entre
// ativacao e calculo de potencial de ignicao e' indiferente, o que elimina
// uma barreira inteira por passo. Os passos 3, 4 e 5 sao feitos por uma
// unica thread (single), pois sao O(1) e mexem em variaveis de controle
// compartilhadas (passo, continua, total, pico).
void simular(Simulacao *sim, int *passos_executados, int *total_ignicoes,
             int *pico_passo, int *pico_qtd) {
    long long n = (long long)sim->L * (long long)sim->C;
    int passo = 0;
    int total = 0;
    int p_passo = -1; // -1 indica "nenhuma ignicao ocorreu" (secao 10)
    int p_qtd = 0;
    // Acumuladores compartilhados usados pelas reducoes. Eles comecam em zero
    // e sao zerados novamente, pela thread do single final, antes do proximo
    // passo. Assim evitamos uma regiao single e uma barreira extras por passo.
    int ignicoes_no_passo = 0;
    int tem_chamas = 0;

    // se nao houver nenhuma celula em chamas logo apos a inicializacao, ou
    // se P for 0, nenhum passo e' executado (secao 9)
    int continua = (sim->P > 0) && existe_em_chamas(sim);
    if (!continua) {
        *passos_executados = 0;
        *total_ignicoes = 0;
        *pico_passo = -1;
        *pico_qtd = 0;
        return;
    }

    #pragma omp parallel num_threads(sim->T) default(none) \
        shared(sim, n, passo, continua, total, p_passo, p_qtd, ignicoes_no_passo, tem_chamas)
    {
        while (continua) {
            // 1. ativar zonas programadas para este passo (inlineado) e
            // 2. calcular proximo estado (cada idx e' independente; so os
            // contadores agregados precisam de reducao, sem critical/atomic)
            // A atualizacao tem custo irregular (celulas intactas examinam
            // oito vizinhos; as demais sao baratas). A distribuicao desse
            // custo muda ao longo da simulacao, a medida que a frente de
            // fogo avanca e as regioes queimadas/contidas crescem, entao um
            // schedule(static) tende a desbalancear conforme T cresce. Entre
            // as tres politicas comparadas experimentalmente (Secao
            // Escalonamento do relatorio: static, dynamic(1024) e
            // guided(1024)), guided(1024) obteve o melhor equilibrio entre
            // balanceamento e overhead de escalonador nas cargas media e
            // grande com T alto, e ficou proximo do melhor tambem na carga
            // pequena; por isso foi escolhida como a politica final. Os
            // blocos comecam grandes e diminuem conforme o trabalho
            // restante; o tamanho minimo de 1024 iteracoes limita a
            // quantidade de requisicoes ao escalonador e preserva a
            // localidade de memoria.
            #pragma omp for schedule(guided, 1024) reduction(+:ignicoes_no_passo) reduction(||:tem_chamas)
            for (long long idx = 0; idx < n; idx++) {
                if (sim->ativacao[idx] == passo && sim->estado_atual[idx] == 1) {
                    sim->estado_atual[idx] = 4;
                }

                int ficou_em_chamas;
                int ignizou = atualizar_celula(sim, idx, &ficou_em_chamas);
                ignicoes_no_passo += ignizou;
                if (ficou_em_chamas) tem_chamas = 1;
            }

            #pragma omp single
            {
                // 3. estatisticas do proximo estado (empate mantem o 1o passo)
                total += ignicoes_no_passo;
                if (ignicoes_no_passo > p_qtd) {
                    p_qtd = ignicoes_no_passo;
                    p_passo = passo;
                }
                trocar_matrizes(sim); // 4. trocar as matrizes
                passo++;
                continua = (passo < sim->P) && tem_chamas; // 5. condicao de parada

                // Prepara as variaveis de reducao para o proximo passo. A
                // barreira implicita do single publica os zeros antes de as
                // threads voltarem ao inicio do while.
                ignicoes_no_passo = 0;
                tem_chamas = 0;
            } // barreira implicita no fim do single: todas as threads veem o novo 'continua'
        }
    }

    *passos_executados = passo;
    *total_ignicoes = total;
    *pico_passo = p_passo;
    *pico_qtd = p_qtd;
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

int main(int argc, char *argv[]) {

    if (argc != 2) {
        printf("Erro de passagem de argumento.\n Uso: ./fire_omp entrada.txt");
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

    int passos_executados, total_ignicoes, pico_passo, pico_qtd;

    // Garante que a regiao solicite exatamente T threads mesmo se o ambiente
    // externo tiver habilitado ajuste dinamico do tamanho da equipe.
    omp_set_dynamic(0);

    double t_inicio = omp_get_wtime();
    simular(&sim, &passos_executados, &total_ignicoes, &pico_passo, &pico_qtd);
    double t_fim = omp_get_wtime();

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