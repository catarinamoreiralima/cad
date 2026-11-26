#!/usr/bin/env bash
# Executa o benchmark comparativo das politicas de escalonamento do laco de
# atualizacao em fire_omp.c (guided(1024), static e dynamic(1024) por
# padrao), chamando scripts/run_macos.sh uma vez para cada uma. Cada
# chamada gera seu proprio CSV em resultados/ (10 repeticoes por
# combinacao de carga e T, T em {2,4,6,8} por padrao -- ver
# scripts/run_macos.sh), exatamente como se o "schedule(...)" da diretiva
# "omp for" principal fosse editado a mao e o script rodado de novo para
# cada politica.
#
# Uso (a partir da raiz do projeto ou de scripts/):
#   ./scripts/run_politicas_macos.sh
# Variaveis opcionais:
#   POLITICAS="static;dynamic,1024;guided,1024"
#       lista separada por ';' com o texto exato a colocar dentro de
#       "schedule(...)" (padrao: static;dynamic, 1024;guided, 1024)
#   RUNNER=scripts/run_macos.sh
#       script chamado para cada politica (troque para scripts/run_linux.sh
#       no cluster Linux)
#   REPETICOES, THREADS, CARGAS, SAIDA
#       repassadas diretamente a $RUNNER (ver seu proprio cabecalho)
#
# fire_omp.c e' editado temporariamente a cada iteracao e SEMPRE restaurado
# ao final (mesmo em erro ou Ctrl-C, via trap), deixando o repositorio no
# schedule original (guided(1024) na versao final do codigo).

set -euo pipefail

RAIZ="$(cd "$(dirname "$0")/.." && pwd)"
cd "$RAIZ"

ARQUIVO_OMP="fire_omp.c"
RUNNER="${RUNNER:-scripts/run_macos.sh}"

# marcador unico da diretiva "omp for" principal (a mesma usada por
# scripts/run_macos.sh para detectar o SCHEDULE); evita mexer em qualquer
# outra linha do arquivo.
MARCADOR='reduction(+:ignicoes_no_passo)'

if ! grep -qF "$MARCADOR" "$ARQUIVO_OMP"; then
    echo "Erro: nao encontrei a diretiva 'omp for' principal em $ARQUIVO_OMP (marcador: $MARCADOR)." >&2
    exit 1
fi

BACKUP="$(mktemp)"
cp "$ARQUIVO_OMP" "$BACKUP"
restaurar() {
    local status=$?
    cp "$BACKUP" "$ARQUIVO_OMP"
    rm -f "$BACKUP"
    echo "-- fire_omp.c restaurado ao schedule original."
    make -f Makefile.macos -B fire_omp >/dev/null 2>&1 || true
    exit "$status"
}
trap restaurar EXIT INT TERM

# Lista de politicas testadas, na ordem em que aparecem no relatorio.
IFS=';' read -r -a POLITICAS_ARR <<< "${POLITICAS:-static;dynamic, 1024;guided, 1024}"

# Arquiva CSVs de rodadas anteriores para que o glob padrao de
# scripts/plot_resultados.py (resultados/tempos_*.csv) nao misture, na
# mesma media, execucoes de configuracoes de build diferentes (por
# exemplo, antes/depois de mudar os flags de otimizacao do Makefile).
if compgen -G "resultados/tempos_*.csv" > /dev/null; then
    DESTINO="resultados/anteriores/$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$DESTINO"
    mv resultados/tempos_*.csv "$DESTINO"/
    echo "-- CSVs anteriores movidos para $DESTINO/"
fi

for politica in "${POLITICAS_ARR[@]}"; do
    echo "=== Politica: schedule($politica) ==="
    sed -i '' "/$MARCADOR/s/schedule([^)]*)/schedule($politica)/" "$ARQUIVO_OMP"
    grep -n 'omp for schedule' "$ARQUIVO_OMP"
    "$RUNNER"
    echo
done

echo "Todas as politicas executadas. CSVs novos em resultados/tempos_*.csv"
