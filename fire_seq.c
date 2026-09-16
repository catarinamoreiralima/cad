#include <stdio.h>
#include <stdlib.h>

typedef struct {
    int linha;
    int coluna;
} Foco;

typedef struct {
    int passo_ativacao;
    int linha_inicial;
    int coluna_inicial;
    int linha_final;
    int coluna_final;
} ZonaContencao;

typedef struct {
    int L, C, P, T, seed, limiar;
    int vento_linha, vento_coluna, intensidade;

    // quantidade de focos iniciais e zonas de contenção
    int F, Z;
    Foco *focos;         
    ZonaContencao *zonas; 
    
    // ponteiros para as matrizes principais da floresta
    int *estado_atual;
    int *proximo_estado;
    int *tempo_atual;
    int *proximo_tempo;
    int *ativacao; 
} Simulacao;

void ler_entrada(const char *arquivo_entrada, Simulacao *sim) {
    FILE *arquivo = fopen(arquivo_entrada, "r");
    if (arquivo == NULL) {
        fprintf(stderr, "Erro: Nao foi possivel abrir o arquivo %s\n", arquivo_entrada);
        exit(EXIT_FAILURE);
    }

    if (fscanf(arquivo, "%d %d %d %d %d %d", 
               &sim->L, &sim->C, &sim->P, &sim->T, &sim->seed, &sim->limiar) != 6) {
        fprintf(stderr, "Erro na leitura da primeira linha.\n");
        exit(EXIT_FAILURE);
    }

    if (sim->L <= 0 || sim->C <= 0 || sim->P < 0 || sim->T <= 0 || sim->limiar <= 0) {
        fprintf(stderr, "Erro: Parametros gerais invalidos.\n");
        exit(EXIT_FAILURE);
    }

    fscanf(arquivo, "%d %d %d", &sim->vento_linha, &sim->vento_coluna, &sim->intensidade);
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

    fscanf(arquivo, "%d %d", &sim->F, &sim->Z);
    if (sim->F < 0 || sim->Z < 0) {
        fprintf(stderr, "Erro: Quantidade de focos ou zonas nao pode ser negativa.\n");
        exit(EXIT_FAILURE);
    }

    sim->focos = (Foco *)malloc(sim->F * sizeof(Foco));
    sim->zonas = (ZonaContencao *)malloc(sim->Z * sizeof(ZonaContencao));

    for (int i = 0; i < sim->F; i++) {
        fscanf(arquivo, "%d %d", &sim->focos[i].linha, &sim->focos[i].coluna);
        
        if (sim->focos[i].linha < 0 || sim->focos[i].linha >= sim->L ||
            sim->focos[i].coluna < 0 || sim->focos[i].coluna >= sim->C) {
            fprintf(stderr, "Erro: Foco fora da matriz.\n");
            exit(EXIT_FAILURE);
        }
        
        for (int j = 0; j < i; j++) {
            if (sim->focos[i].linha == sim->focos[j].linha && 
                sim->focos[i].coluna == sim->focos[j].coluna) {
                fprintf(stderr, "Erro: Foco repetido encontrado.\n");
                exit(EXIT_FAILURE);
            }
        }
    }

    for (int i = 0; i < sim->Z; i++) {
        fscanf(arquivo, "%d %d %d %d %d", 
               &sim->zonas[i].passo_ativacao, 
               &sim->zonas[i].linha_inicial, &sim->zonas[i].coluna_inicial,
               &sim->zonas[i].linha_final, &sim->zonas[i].coluna_final);
        
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

int main(int argc, char *argv[]){
    
    if (argc != 2) {
        printf("Erro de passagem de argumento.\n Uso: ./fire_seq entrada.txt");
        exit(EXIT_FAILURE);
    }
    Simulacao sim;

    ler_entrada(argv[1], &sim);

    // printf("Matriz: %dx%d, Passos: %d, Threads: %d\n", sim.L, sim.C, sim.P, sim.T);
    free(sim.focos);
    free(sim.zonas);

    return 0;
}