#!/usr/bin/env bash
# Executa fire_seq e fire_omp REPETICOES vezes para cada carga e grava os
# tempos em CSV (macOS, clang da Apple + libomp do Homebrew).
#
# Uso (a partir da raiz do projeto ou de scripts/):
#   ./scripts/run_macos.sh
# Variaveis opcionais:
#   REPETICOES=10           execucoes por configuracao
#   THREADS="1 2 4 8"       quantidades de threads da versao OpenMP (padrao: 1 a 8).
#                           Para cada valor, gera uma copia temporaria da entrada
#                           trocando so o T (4o campo da 1a linha).
#                           THREADS=arquivo usa apenas o T do proprio arquivo.
#   CARGAS="pequena media"  subconjunto das cargas
#   SAIDA=arquivo.csv       caminho do CSV
#   SCHEDULE=rotulo         rotulo da coluna schedule (padrao: lido da clausula
#                           schedule do laco de atualizacao em fire_omp.c)

set -euo pipefail

RAIZ="$(cd "$(dirname "$0")/.." && pwd)"
cd "$RAIZ"

NUCLEOS="$(sysctl -n hw.logicalcpu)"
MAKEFILE="Makefile.macos"

REPETICOES="${REPETICOES:-10}"
THREADS="${THREADS:-1 2 3 4 5 6 7 8}"
CARGAS="${CARGAS:-pequena media grande}"
# Politica do laco de atualizacao, extraida do pragma com a reducao de
# ignicoes; "guided, 1024" vira "guided:1024" para nao quebrar o CSV.
if [ -z "${SCHEDULE:-}" ]; then
    SCHEDULE="$(grep 'omp for' fire_omp.c | grep 'reduction(+:ignicoes_no_passo)' |
        sed -n 's/.*schedule(\([^)]*\)).*/\1/p' | tr -d ' ' | tr ',' ':')"
    [ -n "$SCHEDULE" ] || SCHEDULE="desconhecido"
    if [ "$SCHEDULE" = "runtime" ]; then
        SCHEDULE="runtime=$(printf '%s' "${OMP_SCHEDULE:-padrao}" | tr -d ' ' | tr ',' ':')"
    fi
fi
mkdir -p resultados
SAIDA="${SAIDA:-resultados/tempos_macos_${SCHEDULE//:/-}_$(hostname -s)_$(date +%Y%m%d_%H%M%S).csv}"

echo "Compilando com $MAKEFILE..."
make -f "$MAKEFILE" -B fire_seq fire_omp >/dev/null

TMPDIR_EXEC="$(mktemp -d)"
trap 'rm -rf "$TMPDIR_EXEC"' EXIT

# extrai o valor de um campo "nome: valor" da saida do programa
campo() { awk -v c="$1:" '$1 == c { print $2 }' "$2"; }

executar() { # versao carga threads repeticao entrada
    local versao="$1" carga="$2" threads="$3" rep="$4" entrada="$5"
    local out="$TMPDIR_EXEC/out.txt"
    if ! "./fire_$versao" "$entrada" >"$out"; then
        echo "Erro ao executar fire_$versao ($carga, T=$threads)" >&2
        exit 1
    fi
    local sched="$SCHEDULE"
    [ "$versao" = seq ] && sched="-"
    echo "$carga,$versao,$sched,$threads,$rep,$(campo tempo "$out"),$(campo passos "$out"),$(campo checksum "$out")" >>"$SAIDA"
    printf '  %-3s %-7s T=%-3s rep=%s  %ss\n' "$versao" "$carga" "$threads" "$rep" "$(campo tempo "$out")"
}

echo "carga,versao,schedule,threads,repeticao,tempo_s,passos,checksum" >"$SAIDA"
echo "Schedule: $SCHEDULE | Nucleos: $NUCLEOS | Threads: $THREADS | Repeticoes: $REPETICOES"

for carga in $CARGAS; do
    entrada="entrada_carga_$carga.txt"
    [ -f "$entrada" ] || { echo "Arquivo $entrada nao encontrado" >&2; exit 1; }

    for rep in $(seq 1 "$REPETICOES"); do
        executar seq "$carga" 1 "$rep" "$entrada"
    done

    # T e' lido do arquivo (4o campo da 1a linha): gera uma copia por quantidade
    # de threads, ou roda a entrada original com THREADS=arquivo.
    if [ "$THREADS" = arquivo ]; then
        t="$(awk 'NR == 1 { print $4 }' "$entrada")"
        for rep in $(seq 1 "$REPETICOES"); do
            executar omp "$carga" "$t" "$rep" "$entrada"
        done
    else
        for t in $THREADS; do
            entrada_t="$TMPDIR_EXEC/${carga}_T$t.txt"
            awk -v t="$t" 'NR == 1 { $4 = t } { print }' "$entrada" >"$entrada_t"
            for rep in $(seq 1 "$REPETICOES"); do
                executar omp "$carga" "$t" "$rep" "$entrada_t"
            done
        done
    fi
done

# confere se todas as execucoes de uma carga produziram o mesmo checksum
awk -F, 'NR > 1 {
    if (!($1 in ref)) ref[$1] = $8
    else if (ref[$1] != $8) { print "AVISO: checksum divergente em " $1 " (" $2 ", T=" $4 ")"; erro = 1 }
} END { if (!erro) print "Checksums consistentes entre versoes e threads." }' "$SAIDA"

echo "Resultados em $SAIDA"
