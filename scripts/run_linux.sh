#!/usr/bin/env bash
# Executa fire_seq e fire_mpi_omp REPETICOES vezes para cada carga e grava os
# tempos em CSV (Linux, mpicc/mpirun).
#
# Uso (a partir da raiz do projeto ou de scripts/):
#   ./scripts/run_linux.sh
# Variaveis opcionais:
#   REPETICOES=10           execucoes por configuracao
#   PROCESSOS="1 2 4"       quantidades de processos MPI
#   THREADS="1 2 4"         quantidades de threads OpenMP por processo.
#                           Para cada valor, gera uma copia temporaria da entrada
#                           trocando so o T (4o campo da 1a linha).
#                           THREADS=arquivo usa apenas o T do proprio arquivo.
#                           A combinacao P=1 x T=1 da versao MPI e' sempre pulada.
#   MAX_TRABALHADORES=8     limite de P*T; configuracoes acima dele sao puladas
#                           (padrao: nucleos da maquina, no maximo 8)
#   MPIRUN_ARGS="..."       opcoes do mpirun (padrao: "--bind-to none", Open MPI),
#                           ex.: "--bind-to none --host hal01,hal02 --map-by node".
#                           Com varias maquinas, TMPDIR deve apontar para um
#                           diretorio compartilhado, salvo com THREADS=arquivo.
#   CARGAS="pequena media"  subconjunto das cargas
#   SAIDA=arquivo.csv       caminho do CSV
#   SCHEDULE=rotulo         rotulo da coluna schedule (padrao: lido da clausula
#                           schedule do laco de atualizacao em fire_mpi_omp.c)

set -euo pipefail

RAIZ="$(cd "$(dirname "$0")/.." && pwd)"
cd "$RAIZ"

NUCLEOS="$(nproc)"
MAX_TRABALHADORES="${MAX_TRABALHADORES:-$((NUCLEOS < 8 ? NUCLEOS : 8))}"
MAKEFILE="Makefile"

REPETICOES="${REPETICOES:-10}"
PROCESSOS="${PROCESSOS:-1 2 4}"
THREADS="${THREADS:-1 2 4}"
MPIRUN_ARGS="${MPIRUN_ARGS-"--bind-to none"}"
CARGAS="${CARGAS:-pequena media grande}"
# Politica do laco de atualizacao, extraida do pragma com a reducao de
# ignicoes; "guided, 1024" vira "guided:1024" para nao quebrar o CSV.
if [ -z "${SCHEDULE:-}" ]; then
    SCHEDULE="$(grep 'omp for' fire_mpi_omp.c | grep 'reduction(+:ignicoes_no_passo' |
        sed -n 's/.*schedule(\([^)]*\)).*/\1/p' | tr -d ' ' | tr ',' ':')"
    [ -n "$SCHEDULE" ] || SCHEDULE="desconhecido"
    if [ "$SCHEDULE" = "runtime" ]; then
        SCHEDULE="runtime=$(printf '%s' "${OMP_SCHEDULE:-padrao}" | tr -d ' ' | tr ',' ':')"
    fi
fi
mkdir -p resultados
SAIDA="${SAIDA:-resultados/mpi_linux_${SCHEDULE//:/-}_$(hostname -s)_$(date +%Y%m%d_%H%M%S).csv}"

echo "Compilando com $MAKEFILE..."
make -f "$MAKEFILE" -B fire_seq fire_mpi_omp >/dev/null

TMPDIR_EXEC="$(mktemp -d)"
trap 'rm -rf "$TMPDIR_EXEC"' EXIT

# extrai o valor de um campo "nome: valor" da saida do programa
campo() { awk -v c="$1:" '$1 == c { print $2 }' "$2"; }

executar() { # versao carga processos threads repeticao entrada
    local versao="$1" carga="$2" procs="$3" threads="$4" rep="$5" entrada="$6"
    local out="$TMPDIR_EXEC/out.txt"
    local sched="$SCHEDULE"
    local comando=(mpirun -np "$procs" $MPIRUN_ARGS ./fire_mpi_omp "$entrada")
    if [ "$versao" = seq ]; then
        sched="-"
        comando=(./fire_seq "$entrada")
    fi
    if ! "${comando[@]}" >"$out"; then
        echo "Erro ao executar ($carga, P=$procs, T=$threads): ${comando[*]}" >&2
        exit 1
    fi
    echo "$carga,$versao,$sched,$procs,$threads,$rep,$(campo tempo "$out"),$(campo passos "$out"),$(campo checksum "$out")" >>"$SAIDA"
    printf '  %-7s %-7s P=%-3s T=%-3s rep=%s  %ss\n' "$versao" "$carga" "$procs" "$threads" "$rep" "$(campo tempo "$out")"
}

echo "carga,versao,schedule,processos,threads,repeticao,tempo_s,passos,checksum" >"$SAIDA"
echo "Schedule: $SCHEDULE | Nucleos: $NUCLEOS | Limite P*T: $MAX_TRABALHADORES | Processos: $PROCESSOS | Threads: $THREADS | Repeticoes: $REPETICOES"

for carga in $CARGAS; do
    entrada="entrada_carga_$carga.txt"
    [ -f "$entrada" ] || { echo "Arquivo $entrada nao encontrado" >&2; exit 1; }

    for rep in $(seq 1 "$REPETICOES"); do
        executar seq "$carga" 1 1 "$rep" "$entrada"
    done

    # T e' lido do arquivo (4o campo da 1a linha): gera uma copia por quantidade
    # de threads, ou roda a entrada original com THREADS=arquivo.
    lista_t="$THREADS"
    [ "$THREADS" = arquivo ] && lista_t="$(awk 'NR == 1 { print $4 }' "$entrada")"
    for p in $PROCESSOS; do
        for t in $lista_t; do
            # P=1 x T=1 nao e' medido: a referencia de um trabalhador e' fire_seq.
            if [ "$p" -eq 1 ] && [ "$t" -eq 1 ]; then
                continue
            fi
            if [ $((p * t)) -gt "$MAX_TRABALHADORES" ]; then
                echo "Pulando P=$p x T=$t: excede o limite de $MAX_TRABALHADORES"
                continue
            fi
            entrada_t="$entrada"
            if [ "$THREADS" != arquivo ]; then
                entrada_t="$TMPDIR_EXEC/${carga}_T$t.txt"
                awk -v t="$t" 'NR == 1 { $4 = t } { print }' "$entrada" >"$entrada_t"
            fi
            for rep in $(seq 1 "$REPETICOES"); do
                executar mpi_omp "$carga" "$p" "$t" "$rep" "$entrada_t"
            done
        done
    done
done

# confere se todas as execucoes de uma carga produziram o mesmo checksum
awk -F, 'NR > 1 {
    if (!($1 in ref)) ref[$1] = $9
    else if (ref[$1] != $9) { print "AVISO: checksum divergente em " $1 " (" $2 ", P=" $4 ", T=" $5 ")"; erro = 1 }
} END { if (!erro) print "Checksums consistentes entre versoes, processos e threads." }' "$SAIDA"

echo "Resultados em $SAIDA"
