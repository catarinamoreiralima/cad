#!/usr/bin/env python3
"""Gera os graficos de tempo, speedup e eficiencia da versao MPI + OpenMP.

Le os CSVs gerados por scripts/run_linux.sh (colunas: carga,versao,schedule,
processos,threads,repeticao,tempo_s,passos,checksum), usa o tempo medio de cada
configuracao e calcula, com P processos e T threads por processo:

    S(P,T) = t_seq / t_mpi_omp(P,T)
    E(P,T) = S(P,T) / (P*T) * 100%

Produz, para cada schedule e cada quantidade de processos P presente nos dados,
um painel por carga em funcao do numero de threads por processo T:

    relatorio/figuras/tempos_mpi_<schedule>_P<P>.pdf
    relatorio/figuras/speedup_mpi_<schedule>_P<P>.pdf
    relatorio/figuras/eficiencia_mpi_<schedule>_P<P>.pdf

e um resumo tabular em resultados/resumo_mpi.csv. CSVs que casem com o mesmo
padrao sao agregados: nao misture execucoes de ambientes diferentes.

Uso:
    python3 scripts/plot_resultados.py
    python3 scripts/plot_resultados.py --csv "resultados/mpi_linux_*_hal03_*.csv"
    python3 scripts/plot_resultados.py --outdir relatorio/figuras --resumo resultados/resumo_mpi.csv
"""
from __future__ import annotations

import argparse
import csv
import glob
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
MARKERS = ["o", "s", "^", "D", "v", "P"]

CARGA_ORDEM = ["pequena", "media", "grande"]
RAIZ = Path(__file__).resolve().parent.parent


def ler_csvs(padroes: list[str]) -> list[dict]:
    linhas = []
    arquivos = []
    for padrao in padroes:
        arquivos.extend(sorted(glob.glob(padrao)))
    if not arquivos:
        raise SystemExit(f"Nenhum CSV encontrado para: {padroes}")
    for caminho in arquivos:
        with open(caminho, newline="", encoding="utf-8") as f:
            leitor = csv.DictReader(f)
            if "processos" not in (leitor.fieldnames or []):
                raise SystemExit(f"{caminho}: sem a coluna 'processos' (nao e um CSV da versao MPI).")
            for row in leitor:
                row["processos"] = int(row["processos"])
                row["threads"] = int(row["threads"])
                row["tempo_s"] = float(row["tempo_s"])
                linhas.append(row)
    return linhas


def media(valores: list[float]) -> float:
    return sum(valores) / len(valores)


def agrupar(linhas: list[dict]):
    """Agrupa tempo_s por (carga, versao, schedule, processos, threads)."""
    grupos: dict[tuple, list[float]] = defaultdict(list)
    for r in linhas:
        chave = (r["carga"], r["versao"], r["schedule"], r["processos"], r["threads"])
        grupos[chave].append(r["tempo_s"])
    return grupos


def construir_resumo(grupos: dict[tuple, list[float]]):
    """Retorna:
    - seq_base[carga] = tempo medio sequencial
    - pontos = lista de dicts com as metricas de cada configuracao MPI + OpenMP
    """
    seq_base: dict[str, float] = {}
    for (carga, versao, *_resto), tempos in grupos.items():
        if versao == "seq":
            seq_base[carga] = media(tempos)

    pontos = []
    for (carga, versao, schedule, processos, threads), tempos in grupos.items():
        if versao != "mpi_omp" or carga not in seq_base:
            continue
        tempo = media(tempos)
        speedup = seq_base[carga] / tempo
        pontos.append({
            "carga": carga, "schedule": schedule,
            "processos": processos, "threads": threads,
            "total": processos * threads, "tempo": tempo, "speedup": speedup,
            "eficiencia": speedup / (processos * threads) * 100.0, "n": len(tempos),
        })

    pontos.sort(key=lambda p: (ordem_carga(p["carga"]), p["schedule"], p["processos"], p["threads"]))
    return seq_base, pontos


def ordem_carga(carga: str) -> tuple:
    idx = CARGA_ORDEM.index(carga) if carga in CARGA_ORDEM else len(CARGA_ORDEM)
    return (idx, carga)


def escrever_resumo(caminho: Path, seq_base: dict[str, float], pontos: list[dict]) -> None:
    caminho.parent.mkdir(parents=True, exist_ok=True)
    with open(caminho, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["carga", "schedule", "processos", "threads", "total", "t_seq_s",
                    "t_mpi_omp_s", "speedup", "eficiencia_pct", "n_rep"])
        for p in pontos:
            w.writerow([p["carga"], p["schedule"], p["processos"], p["threads"], p["total"],
                        f"{seq_base[p['carga']]:.6f}", f"{p['tempo']:.6f}",
                        f"{p['speedup']:.3f}", f"{p['eficiencia']:.2f}", p["n"]])
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


def plot_metrica(seq_base, pontos, metrica: str, titulo: str, rotulo_y: str, caminho: Path) -> None:
    """P fixo (todos os pontos tem o mesmo): um painel por carga, T no eixo x."""
    cargas = sorted({p["carga"] for p in pontos}, key=ordem_carga)
    n_proc = pontos[0]["processos"]
    threads = sorted({p["threads"] for p in pontos})
    # Tempos tem escalas diferentes entre cargas; speedup e eficiencia sao comparaveis.
    fig, eixos = plt.subplots(1, len(cargas), figsize=(max(7, 4 * len(cargas)), 4.2), dpi=150,
                              sharey=metrica != "tempo", squeeze=False)
    legenda = {}

    for ax, carga in zip(eixos[0], cargas):
        serie = sorted((p for p in pontos if p["carga"] == carga), key=lambda p: p["threads"])
        ax.plot([p["threads"] for p in serie], [p[metrica] for p in serie],
                color=CATEGORICAL[0], marker=MARKERS[0], linewidth=2, markersize=8,
                label="MPI + OpenMP")

        referencia = dict(color=CHROME["muted"], linestyle=":", linewidth=1.5)
        if metrica == "tempo":
            ax.axhline(seq_base[carga], label="Sequencial", **referencia)
            ax.set_ylim(bottom=0, top=max(ax.get_ylim()[1], seq_base[carga] * 1.1))
        elif metrica == "speedup":
            # Amostra a reta S = P*T, que aparece curva com o eixo x logaritmico.
            xs = [threads[0] + (threads[-1] - threads[0]) * k / 100 for k in range(101)]
            ax.plot(xs, [n_proc * x for x in xs], label="Ideal (P x T)", **referencia)
        else:
            ax.axhline(100, label="Ideal", **referencia)

        ax.set_xscale("log", base=2)
        ax.set_xticks(threads, [str(t) for t in threads])
        ax.minorticks_off()
        ax.set_title(f"Carga {carga}")
        ax.set_xlabel("Threads por processo (T)")
        estilizar_eixo(ax)
        for linha, rotulo in zip(*ax.get_legend_handles_labels()):
            legenda.setdefault(rotulo, linha)

    eixos[0][0].set_ylabel(rotulo_y)
    if metrica != "tempo":
        eixos[0][0].set_ylim(bottom=0)
    fig.suptitle(titulo, color=CHROME["text_primary"])
    fig.legend(legenda.values(), legenda.keys(), frameon=False, fontsize=9,
               loc="lower center", ncol=min(len(legenda), 5))
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(caminho, facecolor=fig.get_facecolor())
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv", nargs="+", default=[str(RAIZ / "resultados" / "mpi_*.csv")],
                         help="Padrao(oes) glob para os CSVs de entrada (padrao: resultados/mpi_*.csv)")
    parser.add_argument("--outdir", default=str(RAIZ / "relatorio" / "figuras"),
                         help="Diretorio de saida para os PDFs (padrao: relatorio/figuras)")
    parser.add_argument("--resumo", default=str(RAIZ / "resultados" / "resumo_mpi.csv"),
                         help="Caminho do CSV-resumo (padrao: resultados/resumo_mpi.csv)")
    args = parser.parse_args()

    linhas = ler_csvs(args.csv)
    seq_base, pontos = construir_resumo(agrupar(linhas))

    if not seq_base:
        raise SystemExit("Nenhuma execucao 'seq' encontrada nos CSVs; nao e possivel calcular speedup/eficiencia.")
    if not pontos:
        raise SystemExit("Nenhuma execucao 'mpi_omp' encontrada nos CSVs.")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    metricas = [
        ("tempos", "tempo", "Tempo de execucao", "Tempo (s)"),
        ("speedup", "speedup", "Speedup", "Speedup S(P,T)"),
        ("eficiencia", "eficiencia", "Eficiencia", "Eficiencia E(P,T) (%)"),
    ]
    for schedule, n_proc in sorted({(p["schedule"], p["processos"]) for p in pontos}):
        selecao = [p for p in pontos if p["schedule"] == schedule and p["processos"] == n_proc]
        sufixo = "_mpi_" + schedule.replace(":", "-").replace("=", "-") + f"_P{n_proc}"
        for prefixo, metrica, titulo, rotulo_y in metricas:
            caminho = outdir / f"{prefixo}{sufixo}.pdf"
            plot_metrica(seq_base, selecao, metrica,
                         f"{titulo} - {n_proc} processo(s) MPI ({schedule})", rotulo_y, caminho)
            print(f"Grafico escrito em {caminho}")

    escrever_resumo(Path(args.resumo), seq_base, pontos)


if __name__ == "__main__":
    main()
