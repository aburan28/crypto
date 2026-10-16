from pathlib import Path
import csv

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


ROOT = Path(__file__).resolve().parent
rows = list(csv.DictReader((ROOT / "data.tsv").open(), delimiter="\t"))
labels = [r["variant"] for r in rows]
median = [float(r["median_bps"]) for r in rows]
low = [float(r["min_bps"]) for r in rows]
high = [float(r["max_bps"]) for r in rows]
colors = ["#4472c4", "#a5a5a5", "#70ad47"]


def chart(ax):
    x = list(range(len(rows)))
    yerr = [[m - lo for m, lo in zip(median, low)], [hi - m for hi, m in zip(high, median)]]
    ax.bar(x, median, yerr=yerr, capsize=5, color=colors, edgecolor="#333333", linewidth=0.7)
    ax.axhline(20.0, color="#c00000", linestyle="--", linewidth=1.4, label="20 B/s target")
    for i, value in enumerate(median):
        ax.text(i, value + 0.45, f"{value:.3f}", ha="center", va="bottom", fontsize=10)
    ax.set_xticks(x, ["live\n256 threads", "matched control\n512 threads", "candidate\n512 + table"])
    ax.set_ylim(0, 21.8)
    ax.set_ylabel("B/s (billion iterations per second)")
    ax.set_title("DP32 ordinary sigma walk — RunPod RTX PRO 6000")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(loc="upper right", frameon=False)


fig, ax = plt.subplots(figsize=(10.5, 5.8), constrained_layout=True)
chart(ax)
fig.text(0.01, 0.01, "Median of n=3 rotating-panel samples; whiskers show min/max. Run 20261012T1815Z.", fontsize=8)
fig.savefig(ROOT / "dp32-rates.svg", format="svg")
plt.close(fig)


with PdfPages(ROOT / "report.pdf") as pdf:
    page = plt.figure(figsize=(11.7, 8.3))
    page.text(0.06, 0.92, "DP32 ordinary-sigma table screen", fontsize=22, weight="bold")
    page.text(0.06, 0.875, "RunPod · NVIDIA RTX PRO 6000 Blackwell Server Edition · 20261012T1815Z", fontsize=11)
    page.text(0.06, 0.81, "Decision", fontsize=15, weight="bold")
    page.text(
        0.06,
        0.765,
        "The exact polynomial-coordinate sigma table passed the DP32 correctness screen,\n"
        "but did not reach 20 B/s or beat the live geometry on this card. Candidate: 9.254351 B/s\n"
        "(+0.7947% vs matched 512-thread control; -1.0330% vs live 256-thread geometry).",
        fontsize=12,
        va="top",
    )
    table_data = [
        ["Arm", "Median B/s", "Range B/s", "Table"],
        ["live / 256 threads", "9.350950", "9.348547–9.354615", "off"],
        ["control / 512 threads", "9.181391", "9.179768–9.181455", "off"],
        ["candidate / 512 threads", "9.254351", "9.252145–9.255851", "on"],
    ]
    table = page.add_axes([0.06, 0.48, 0.88, 0.19])
    table.axis("off")
    t = table.table(cellText=table_data[1:], colLabels=table_data[0], loc="center", cellLoc="left")
    t.auto_set_font_size(False)
    t.set_fontsize(11)
    t.scale(1, 1.65)
    page.text(0.06, 0.40, "Correctness receipt", fontsize=15, weight="bold")
    page.text(
        0.06,
        0.355,
        "• CUDA arithmetic suite passed for all arms\n"
        "• 5/5 DP32 points reference-verified per arm; 0 dropped\n"
        "• sorted v2 corpora identical across live, control, and candidate\n"
        "• candidate: 126 registers/thread, 68,352 dynamic shared bytes, one block/SM\n"
        "• this is a screen, not full campaign qualification; production was unchanged",
        fontsize=11,
        va="top",
    )
    page.text(0.06, 0.12, "Source: source-files.sha256 and raw/run11; prior 300-report timeout is in raw/failed-run9.", fontsize=9)
    pdf.savefig(page, bbox_inches="tight")
    plt.close(page)

    graph_page, graph_ax = plt.subplots(figsize=(11.7, 8.3))
    chart(graph_ax)
    graph_page.text(0.06, 0.02, "Whiskers show the three recorded A/B samples; dashed line is the requested 20 B/s target.", fontsize=9)
    pdf.savefig(graph_page, bbox_inches="tight")
    plt.close(graph_page)
