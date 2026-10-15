"""Lay out the recorded values of sessions A-D (no computation beyond means
of recorded phase costs).  Run with a Python that has matplotlib:

    python research/ecbench_rootfind_speedup_20261014/figures/plot.py
"""

import json
import pathlib
import statistics

import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
STUDY = HERE.parent
plt.rcParams["svg.hashsalt"] = "ecbench-rootfind-speedup"
plt.rcParams["font.size"] = 9

SESSIONS = {
    "A": "claim-n13",
    "B": "check-n19",
    "C": "claim-renew-n13",
    "D": "claim-renew-n19",
}
COLOURS = {
    "walk": "#7f7f7f",
    "stride": "#4c78a8",
    "renew": "#2a9d4b",
    "rho": "#d1495b",
}


def table_rows(name):
    """Arm, mean S, interval and lower-bound flag as `ecbench table` printed them."""
    rows = {}
    for line in (STUDY / "tables" / f"table-{name}.md").read_text().splitlines():
        cells = [c.strip() for c in line.split("|")]
        if len(cells) < 15 or cells[1] in ("session", "---") or cells[1].startswith("-"):
            continue
        lo, hi = cells[8].strip("[]").split(",")
        rows[cells[2]] = (float(cells[7]), float(lo), float(hi), cells[13])
    return rows


def session_id(name):
    return json.loads((STUDY / "sessions" / name / "session.json").read_text())["session_id"]


def colour(arm):
    for key in ("walk", "stride", "renew"):
        if arm.endswith(key):
            return COLOURS[key]
    return COLOURS["rho"]


def s_figure(panels=None, out=HERE / "s-by-arm"):
    panels = panels or (
        ("n = 13 (r = 2003): sessions A and C", ("A", "C")),
        ("n = 19 (r = 130873): sessions B and D", ("B", "D")),
    )
    fig, axes = plt.subplots(1, len(panels), figsize=(10, 4.2), sharey=True)
    for ax, (title, names) in zip(axes, panels):
        labels, xs = [], 0
        for key in names:
            rows = table_rows(SESSIONS[key])
            for arm, (mean, lo, hi, lb) in rows.items():
                if arm.endswith("-aa"):
                    continue
                ax.errorbar(
                    xs,
                    mean,
                    yerr=[[mean - lo], [hi - mean]],
                    fmt="o" if lb == "no" else "^",
                    color=colour(arm),
                    capsize=3,
                )
                labels.append(f"{key}: {arm.replace('ic-negfold-m2-', '')}")
                xs += 1
            xs += 0.5
        ticks = [i for i, _ in enumerate(labels)]
        pos, k = [], 0
        for key in names:
            n = len([a for a in table_rows(SESSIONS[key]) if not a.endswith("-aa")])
            pos += [k + i for i in range(n)]
            k += n + 0.5
        ax.set_xticks(pos, labels, rotation=55, ha="right")
        ax.set_yscale("log")
        ax.set_title(title)
        ax.grid(axis="y", alpha=0.3)
    axes[0].set_ylabel("mean S = GAE / sqrt(r)  (log scale; bars: 95% interval)")
    fig.suptitle(
        "Counted cost of complete IC (m = 2, rootfind-negfold) by target source, against rho\n"
        "circle = fully priced; triangle = lower bound (unpriced work); "
        + ", ".join(f"{k} {session_id(SESSIONS[k])}" for _, ks in panels for k in ks),
        fontsize=8,
    )
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out.with_suffix(".svg"))
    fig.savefig(out.with_suffix(".png"), dpi=150)
    plt.close(fig)


def phase_means(name, arm):
    out = {}
    for line in (STUDY / "sessions" / name / "records.jsonl").read_text().splitlines():
        rec = json.loads(line)
        if rec["arm"] != arm or rec["warmup"]:
            continue
        for phase in rec["phases"]:
            out.setdefault(phase["name"], []).append(phase["gae"])
    return {k: statistics.mean(v) for k, v in out.items()}


def phase_figure():
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8))
    phases = ("oracle_setup", "relations", "verify", "factor_base", "linear_algebra")
    shades = ("#b07aa1", "#4c78a8", "#f2a541", "#9c9c9c", "#59a14f")
    for ax, (key, arms) in zip(
        axes,
        (
            ("C", ("ic-negfold-m2-d5-walk", "ic-negfold-m2-d5-stride", "ic-negfold-m2-d5-renew", "ic-negfold-m2-d4-renew")),
            ("D", ("ic-negfold-m2-d6-walk", "ic-negfold-m2-d6-stride", "ic-negfold-m2-d6-renew")),
        ),
    ):
        name = SESSIONS[key]
        for i, arm in enumerate(arms):
            means, bottom = phase_means(name, arm), 0.0
            for phase, shade in zip(phases, shades):
                v = means.get(phase, 0.0)
                ax.bar(i, v, bottom=bottom, color=shade, label=phase if i == 0 else None)
                bottom += v
        ax.set_xticks(range(len(arms)), [a.replace("ic-negfold-m2-", "") for a in arms])
        ax.set_ylabel("mean GAE per complete solve")
        ax.set_title(f"session {key} ({session_id(name)})", fontsize=9)
        ax.grid(axis="y", alpha=0.3)
    axes[0].legend(fontsize=7)
    fig.suptitle("Where the counted work goes, by phase (mean of measured runs, warm-ups excluded)", fontsize=9)
    fig.tight_layout()
    fig.savefig(HERE / "phases.svg")
    fig.savefig(HERE / "phases.png", dpi=150)


if __name__ == "__main__":
    s_figure()
    phase_figure()
    # One graph per run report (sessions launched together share a run).
    s_figure(
        (("A: n = 13, stride", ("A",)), ("B: n = 19, stride", ("B",))),
        STUDY / "reports" / "run-20261014T1700Z" / "s-by-arm",
    )
    s_figure(
        (("C: n = 13, stride-renew", ("C",)), ("D: n = 19, stride-renew", ("D",))),
        STUDY / "reports" / "run-20261014T1728Z" / "s-by-arm",
    )
