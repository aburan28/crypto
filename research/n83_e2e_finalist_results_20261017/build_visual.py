#!/usr/bin/env python3
"""Plot the provider-bound sixteen-point Tier 1 finalist result."""

import hashlib
import json
import os
import statistics
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", tempfile.gettempdir())
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

HERE = Path(__file__).resolve().parent
FINAL = HERE.parent / "n83_e2e_final_holdout_20261016"
LABELS = ("1104/r11", "1104/r22", "1107/r11", "1107/r22",
          "1112/r11", "1112/r22", "1120/r11", "1120/r22")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def quartiles(values):
    ordered = sorted(values)
    return [statistics.median(ordered[:8]), statistics.median(ordered),
            statistics.median(ordered[8:])]


def main():
    validated_path = FINAL / "validated_rows.json"
    selection_path = FINAL / "selection.json"
    validated = json.loads(validated_path.read_text())
    selection = json.loads(selection_path.read_text())
    if (validated["evidence_status"] != "PASS_INDEPENDENT_REPLAY"
            or selection["admission_status"] != "PASS"
            or selection["point_count"] != 16
            or len(selection["candidate_ids"]) != len(LABELS)):
        raise ValueError("finalist replay and Tier 1 admission are required")
    ordered = sorted(validated["points"], key=lambda point: point["point_index"])
    if [point["point_index"] for point in ordered] != list(range(16)):
        raise ValueError("sixteen exact public points are required")
    data = []
    for label, candidate in zip(LABELS, selection["candidate_ids"]):
        arms = [next(arm for arm in point["arms"]
                     if arm["candidate_id"] == candidate) for point in ordered]
        online = [arm["one_target_online_ms"] for arm in arms]
        cold = [arm["cold_process_ms"] / 1000 for arm in arms]
        rho_online = [point["rho_by_candidate"][candidate]["online_ms"]
                      for point in ordered]
        ratios = selection["rho_over_ic_online_by_candidate"][candidate]
        if (len(online) != 16 or len(cold) != 16 or len(ratios) != 16
                or len(rho_online) != 16
                or any(value <= 0 for value in online + cold + ratios + rho_online)):
            raise ValueError(f"missing paired finalist measurements: {label}")
        data.append({"label": label, "candidate_id": candidate,
                     "actual_usable_points": arms[0]["raw_replay"]["actual_factor_base_points"],
                     "effective_columns": arms[0]["raw_replay"]["effective_columns"],
                     "online_ms": online, "cold_s": cold,
                     "rho_over_ic_online": ratios,
                     "ic_faster_points": sum(value > 1 for value in ratios),
                     "online_quartiles_ms": quartiles(online),
                     "cold_mean_s": statistics.mean(cold),
                     "rho_mean_online_ms": statistics.mean(rho_online),
                     "cold_mean_over_rho_mean_online":
                     1000 * statistics.mean(cold) / statistics.mean(rho_online),
                     "median_paired_ic_over_rho_online":
                     statistics.median(1 / value for value in ratios),
                     "cold_quartiles_s": quartiles(cold),
                     "ratio_quartiles": quartiles(ratios)})
    record = {"schema": "n83.finalist-plot-data/v1",
              "validated_rows_sha256": digest(validated_path),
              "selection_sha256": digest(selection_path),
              "run_id": validated["run_id"],
              "points": 16, "candidates": data}
    (HERE / "plot_data.json").write_text(json.dumps(record, sort_keys=True,
                                                       indent=2) + "\n")

    plt.rcParams.update({"font.size": 10, "svg.hashsalt": "n83-finalist-39005111289"})
    fig, axes = plt.subplots(1, 3, figsize=(13.4, 5.1))
    fig.subplots_adjust(left=0.055, right=0.99, top=0.82, bottom=0.23,
                        wspace=0.18)
    columns = (("online_ms", "Online IC, milliseconds", True),
               ("cold_s", "Cold IC, seconds", False),
               ("rho_over_ic_online", "Same-point rho / IC online", True))
    for ax, (key, title, logarithmic) in zip(axes, columns):
        for index, item in enumerate(data):
            values = item[key]
            color = "#17688f" if index % 2 == 0 else "#cf752b"
            left, middle, right = quartiles(values)
            ax.vlines(index, left, right, color=color, lw=6, alpha=0.45)
            ax.scatter([index] * len(values), values, s=8, color=color,
                       alpha=0.36, zorder=3)
            ax.scatter(index, middle, marker="_", s=170, linewidths=2.2,
                       color="#172d40", zorder=4)
        ax.set_title(title, fontsize=11, weight="bold")
        ax.set_xticks(range(len(data)), [item["label"] for item in data],
                      rotation=52, ha="right")
        ax.set_xlim(-0.55, 7.55)
        ax.grid(axis="y", alpha=0.25)
        ax.set_axisbelow(True)
        if logarithmic:
            ax.set_yscale("log")
            ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    fig.suptitle("N83 packed-root factor bases · 16 fresh points · two admitted Tier 1 hosts",
                 fontsize=13, weight="bold", y=0.97)
    fig.text(0.5, 0.025,
             "Each dot is one verified public point; thick bar is IQR; dark tick is median. Blue: root11, orange: root22.",
             ha="center", fontsize=9)
    fig.savefig(HERE / "finalist.svg", metadata={"Date": None})
    fig.savefig(HERE / "finalist.png", dpi=170, metadata={"Software": None})
    plt.close(fig)
    print(json.dumps({"status": "PASS", "run_id": validated["run_id"],
                      "candidate_count": len(data)}, sort_keys=True))


if __name__ == "__main__":
    main()
