#!/usr/bin/env python3
"""Gera os graficos de tempo, speedup e eficiencia pedidos no Trab01-Fire-OMP.

Le um ou mais CSVs de resultados/ (colunas: carga,versao,schedule,threads,
repeticao,tempo_s,passos,checksum), calcula as metricas definidas no
relatorio (S(T) = t_seq / t_omp(T), E(T) = S(T)/T * 100%) e produz:

    relatorio/figuras/tempos.pdf
    relatorio/figuras/speedup.pdf
    relatorio/figuras/eficiencia.pdf

e um resumo tabular em resultados/resumo.csv, pronto para preencher a
tabela "Tempos e metricas de desempenho" do relatorio.

Uso:
    python3 scripts/plot_resultados.py
    python3 scripts/plot_resultados.py --csv "resultados/tempos_*.csv"
    python3 scripts/plot_resultados.py --outdir relatorio/figuras --resumo resultados/resumo.csv
"""
from __future__ import annotations

import argparse
import csv
import glob
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Paleta categorica validada (dataviz skill), modo claro.
CHROME = {
    "surface": "#fcfcfb",
    "text_primary": "#0b0b0b",
    "text_secondary": "#52514e",
    "muted": "#898781",
    "grid": "#e1e0d9",
    "axis": "#c3c2b7",
}
CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
LINESTYLES = ["-", "--", "-.", ":"]

CARGA_ORDEM = ["pequena", "media", "grande"]


def carga_cor(carga: str) -> str:
    if carga in CARGA_ORDEM:
        return CATEGORICAL[CARGA_ORDEM.index(carga) % len(CATEGORICAL)]
    # carga desconhecida: cai nos slots seguintes, em ordem alfabetica estavel
    extras = sorted(c for c in {carga} - set(CARGA_ORDEM))
    idx = len(CARGA_ORDEM) + extras.index(carga)
    return CATEGORICAL[idx % len(CATEGORICAL)]


def carga_label(carga: str) -> str:
    return carga[:1].upper() + carga[1:]


def ler_csvs(padroes: list[str]) -> list[dict]:
    linhas = []
    arquivos = []
    for padrao in padroes:
        arquivos.extend(sorted(glob.glob(padrao)))
    if not arquivos:
        raise SystemExit(f"Nenhum CSV encontrado para: {padroes}")
    for caminho in arquivos:
        with open(caminho, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                row["threads"] = int(row["threads"])
                row["tempo_s"] = float(row["tempo_s"])
                linhas.append(row)
    return linhas


def media(valores: list[float]) -> float:
    return sum(valores) / len(valores)


def agrupar(linhas: list[dict]):
    """Agrupa tempo_s por (carga, versao, schedule, threads)."""
    grupos: dict[tuple, list[float]] = defaultdict(list)
    for r in linhas:
        chave = (r["carga"], r["versao"], r["schedule"], r["threads"])
        grupos[chave].append(r["tempo_s"])
    return grupos


def construir_resumo(grupos: dict[tuple, list[float]]):
    """Retorna:
    - seq_base[carga] = tempo medio sequencial
    - serie[(carga, schedule)] = lista de (threads, t_omp_medio, speedup, eficiencia, n)
    """
    seq_base: dict[str, float] = {}
    for (carga, versao, _schedule, _t), tempos in grupos.items():
        if versao == "seq":
            seq_base[carga] = media(tempos)

    serie: dict[tuple[str, str], list[tuple[int, float, float, float, int]]] = defaultdict(list)
    for (carga, versao, schedule, threads), tempos in grupos.items():
        if versao != "omp":
            continue
        if carga not in seq_base:
            continue
        t_omp = media(tempos)
        t_seq = seq_base[carga]
        speedup = t_seq / t_omp
        eficiencia = speedup / threads * 100.0
        serie[(carga, schedule)].append((threads, t_omp, speedup, eficiencia, len(tempos)))

    for chave in serie:
        serie[chave].sort(key=lambda x: x[0])

    return seq_base, serie


def escrever_resumo(caminho: Path, seq_base: dict[str, float], serie: dict) -> None:
    caminho.parent.mkdir(parents=True, exist_ok=True)
    with open(caminho, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["carga", "schedule", "threads", "t_seq_s", "t_omp_s", "speedup", "eficiencia_pct", "n_rep"])
        for (carga, schedule), pontos in sorted(serie.items(), key=ordem_serie):
            for threads, t_omp, speedup, eficiencia, n in pontos:
                w.writerow([carga, schedule, threads, f"{seq_base[carga]:.6f}", f"{t_omp:.6f}",
                            f"{speedup:.3f}", f"{eficiencia:.2f}", n])
    print(f"Resumo escrito em {caminho}")


def estilizar_eixo(ax) -> None:
    ax.set_facecolor(CHROME["surface"])
    ax.figure.set_facecolor(CHROME["surface"])
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(CHROME["axis"])
    ax.tick_params(colors=CHROME["text_secondary"])
    ax.grid(True, color=CHROME["grid"], linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.xaxis.label.set_color(CHROME["text_primary"])
    ax.yaxis.label.set_color(CHROME["text_primary"])
    ax.title.set_color(CHROME["text_primary"])


def linestyle_por_schedule(schedules: list[str]) -> dict[str, str]:
    mapa = {}
    outros = [s for s in schedules if s != "-"]
    for i, s in enumerate(sorted(set(outros))):
        mapa[s] = LINESTYLES[i % len(LINESTYLES)]
    mapa["-"] = ":"
    return mapa


def ordem_serie(item) -> tuple:
    (carga, schedule), _pontos = item
    idx = CARGA_ORDEM.index(carga) if carga in CARGA_ORDEM else 99
    return (idx, carga, schedule)


def plot_tempos(seq_base, serie, outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=150)
    schedules = [s for (_c, s) in serie.keys()]
    ls_map = linestyle_por_schedule(schedules)

    todos_threads = sorted({t for pontos in serie.values() for (t, *_r) in pontos})
    max_t = max(todos_threads) if todos_threads else 8

    for carga in sorted(seq_base, key=lambda c: CARGA_ORDEM.index(c) if c in CARGA_ORDEM else 99):
        cor = carga_cor(carga)
        ax.hlines(seq_base[carga], 1, max_t, colors=cor, linestyles=":", linewidth=1.5, alpha=0.6,
                   label=f"{carga_label(carga)} - sequencial")

    for (carga, schedule), pontos in sorted(serie.items(), key=ordem_serie):
        cor = carga_cor(carga)
        xs = [p[0] for p in pontos]
        ys = [p[1] for p in pontos]
        rotulo_sched = f" ({schedule})" if schedule != "-" else ""
        ax.plot(xs, ys, color=cor, linestyle=ls_map[schedule], linewidth=2, marker="o", markersize=6,
                label=f"{carga_label(carga)} - OpenMP{rotulo_sched}")

    ax.set_xlabel("Threads (T)")
    ax.set_ylabel("Tempo (s)")
    ax.set_title("Tempo de execucao em funcao do numero de threads")
    ax.set_xticks(todos_threads)
    estilizar_eixo(ax)
    ax.legend(frameon=False, fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(outdir / "tempos.pdf", facecolor=fig.get_facecolor())
    plt.close(fig)


def plot_speedup(serie, outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=150)
    schedules = [s for (_c, s) in serie.keys()]
    ls_map = linestyle_por_schedule(schedules)

    todos_threads = sorted({t for pontos in serie.values() for (t, *_r) in pontos})
    max_t = max(todos_threads) if todos_threads else 8

    ax.plot([1, max_t], [1, max_t], color=CHROME["muted"], linestyle="--", linewidth=1.5,
            label="Speedup ideal S(T)=T")

    for (carga, schedule), pontos in sorted(serie.items(), key=ordem_serie):
        cor = carga_cor(carga)
        xs = [p[0] for p in pontos]
        ys = [p[2] for p in pontos]
        rotulo_sched = f" ({schedule})" if schedule != "-" else ""
        ax.plot(xs, ys, color=cor, linestyle=ls_map[schedule], linewidth=2, marker="o", markersize=6,
                label=f"{carga_label(carga)}{rotulo_sched}")

    ax.set_xlabel("Threads (T)")
    ax.set_ylabel("Speedup S(T)")
    ax.set_title("Speedup em funcao do numero de threads")
    ax.set_xticks(todos_threads)
    estilizar_eixo(ax)
    ax.legend(frameon=False, fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(outdir / "speedup.pdf", facecolor=fig.get_facecolor())
    plt.close(fig)


def plot_eficiencia(serie, outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=150)
    schedules = [s for (_c, s) in serie.keys()]
    ls_map = linestyle_por_schedule(schedules)

    todos_threads = sorted({t for pontos in serie.values() for (t, *_r) in pontos})

    ax.axhline(100, color=CHROME["muted"], linestyle="--", linewidth=1.5, label="Eficiencia ideal (100%)")

    for (carga, schedule), pontos in sorted(serie.items(), key=ordem_serie):
        cor = carga_cor(carga)
        xs = [p[0] for p in pontos]
        ys = [p[3] for p in pontos]
        rotulo_sched = f" ({schedule})" if schedule != "-" else ""
        ax.plot(xs, ys, color=cor, linestyle=ls_map[schedule], linewidth=2, marker="o", markersize=6,
                label=f"{carga_label(carga)}{rotulo_sched}")

    ax.set_xlabel("Threads (T)")
    ax.set_ylabel("Eficiencia E(T) (%)")
    ax.set_title("Eficiencia em funcao do numero de threads")
    ax.set_xticks(todos_threads)
    estilizar_eixo(ax)
    ax.legend(frameon=False, fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(outdir / "eficiencia.pdf", facecolor=fig.get_facecolor())
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv", nargs="+", default=["/Users/catarinamoreiralima/Documents/Faculdade/cad/trabalho-01/resultados/tempos_*.csv"],
                         help="Padrao(oes) glob para os CSVs de entrada (padrao: resultados/tempos_*.csv)")
    parser.add_argument("--outdir", default="relatorio/figuras",
                         help="Diretorio de saida para os PDFs (padrao: relatorio/figuras)")
    parser.add_argument("--resumo", default="resultados/resumo.csv",
                         help="Caminho do CSV-resumo com as metricas (padrao: resultados/resumo.csv)")
    args = parser.parse_args()

    linhas = ler_csvs(args.csv)
    grupos = agrupar(linhas)
    seq_base, serie = construir_resumo(grupos)

    if not seq_base:
        raise SystemExit("Nenhuma execucao 'seq' encontrada nos CSVs; nao e possivel calcular speedup/eficiencia.")
    if not serie:
        raise SystemExit("Nenhuma execucao 'omp' encontrada nos CSVs.")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    plot_tempos(seq_base, serie, outdir)
    plot_speedup(serie, outdir)
    plot_eficiencia(serie, outdir)
    escrever_resumo(Path(args.resumo), seq_base, serie)

    print(f"Graficos escritos em {outdir}/tempos.pdf, {outdir}/speedup.pdf, {outdir}/eficiencia.pdf")


if __name__ == "__main__":
    main()
