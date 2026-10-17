"""The repository's standard isogeny-volcano heat map (`volcano-diagram/v2`).

Draws a recorded ℓ-isogeny volcano as stacked strata: each level is a
shaded band holding an elliptical ring seen in perspective, the crater on
top, every lower ring wider than the one above, so the figure reads as a
volcano. Each vertex sits on its ring at the angle of the curves it leads
down to, under its parent, and edges join parents to children and crater
neighbours. One panel per metric (e.g. mean operations and mean time), each
on one observed-range YlOrRd scale shared by every level. Vertices without
a value are small and grey, so coverage is visible. Plotting only: every
number it draws comes from the input (AGENTS.md "Implementation language").

usage:
    python volcano_diagram.py INPUT.json OUT_STEM     # writes OUT_STEM.svg and OUT_STEM.png (150 dpi)

or from Python:
    from volcano_diagram import render
    render(spec, "out/volcano")

Input (`volcano-diagram/v2`, see README.md):
    {
      "schema": "volcano-diagram/v2",
      "ell": 2, "depth": 3,
      "title": "Measured 2-isogeny volcano heat map · ...",
      "subtitle": "run-...",                       # run or session id
      "footnote": "...",                            # optional line under the panels
      "metrics": [
        {"key": "ops",  "title": "Native rho · mean operations", "label": "group operations per solve"},
        {"key": "time", "title": "Native rho · mean time",       "label": "ms per solve"}
      ],
      "vertices": [
        {"id": 0, "label": "v0", "level": 3, "nbrs": [1],
         "values": {"ops": 658723.0, "time": 431.2},   # optional, per metric key
         "highlight": "root"}                          # optional: ringed blue, named in the title
      ],
      "edges": [{"a": 0, "b": 1, "status": "unverified"}]   # optional; default: nbrs, verified
    }
`volcano-diagram/v1` input (a single `value` per vertex and `value_label`)
is still accepted and drawn as one panel. Levels count down from the
crater: 0 = crater, depth = floor. Vertices with `level: null` are left
out. Unverified edges are dashed.
"""

import json
import math
import sys
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Ellipse, Rectangle  # noqa: E402

SCHEMA = "volcano-diagram/v2"
CMAP = "YlOrRd"
EDGE = "#c3c7cf"
BAND = "#f4f6fb"
INK = "#2b3445"
MUTED = "#5f6b7d"
RING = "#1f4e9e"


def _v2(spec):
    """Accept v1 input as a one-metric v2 spec."""
    if spec.get("schema") == SCHEMA:
        return spec
    if spec.get("schema") != "volcano-diagram/v1":
        raise ValueError(f"expected schema {SCHEMA} or volcano-diagram/v1")
    out = dict(spec)
    out["metrics"] = [{"key": "value", "title": spec.get("value_label", "value"), "label": spec.get("value_label", "value")}]
    out["vertices"] = [dict(v, values={"value": v.get("value")}) for v in spec["vertices"]]
    return out


def _crater_order(verts):
    """Crater vertices in cycle order: start at an end of the drawn path (a
    vertex with one crater neighbour) or, for a full cycle, the least id,
    and follow horizontal edges; remaining pieces follow in id order."""
    crater = {v for v, d in verts.items() if d["level"] == 0}
    hn = {v: sorted(w for w in verts[v]["nbrs"] if w in crater and w != v) for v in crater}
    order, seen = [], set()
    while len(seen) < len(crater):
        rest = sorted(crater - seen)
        ends = [v for v in rest if len([w for w in hn[v] if w not in seen]) <= 1]
        v = ends[0] if ends else rest[0]
        while v is not None:
            order.append(v)
            seen.add(v)
            v = next((w for w in hn[v] if w not in seen), None)
    return order


def _positions(verts):
    """Position around the volcano, in [0, 1): leaves of each crater
    vertex's tree consecutive, parents at their children's mean, so every
    vertex sits at the angle of the curves below it."""
    crater = _crater_order(verts)
    children = defaultdict(list)
    for v, d in sorted(verts.items()):
        if d["level"] > 0:
            ups = sorted(w for w in d["nbrs"] if w in verts and verts[w]["level"] == d["level"] - 1)
            if ups:
                children[ups[0]].append(v)
    pos, cursor = {}, [0.0]

    def place(v):
        if not children[v]:
            pos[v] = cursor[0]
            cursor[0] += 1.0
        else:
            for c in children[v]:
                place(c)
            pos[v] = sum(pos[c] for c in children[v]) / len(children[v])

    for c in crater:
        place(c)
        cursor[0] += 0.5
    # Vertices whose tree never reaches the crater (a truncated walk) go last.
    for v in sorted(verts, key=lambda v: (verts[v]["level"], v)):
        if v not in pos:
            place(v)
    total = cursor[0] or 1.0
    return {v: p / total for v, p in pos.items()}, children


def _geometry(depth):
    """Ring radii and centres per level: wider and lower going down."""
    rx = [1.0 + 0.32 * i for i in range(depth + 1)]
    ry = [0.30 * r for r in rx]
    yc = [0.0]
    for i in range(1, depth + 1):
        yc.append(yc[-1] - (ry[i - 1] + ry[i] + 0.55))
    return rx, ry, yc


def _fmt(x):
    a = abs(x)
    if a >= 1e5:
        return f"{x:,.0f}"
    if a >= 100:
        return f"{x:,.2f}"
    return f"{x:.3f}"


def _panel(ax, spec, metric, verts, pos, edges, depth):
    rx, ry, yc = _geometry(depth)
    width = rx[-1] + 0.55
    # Turn the volcano so the highlighted vertex (the walk's root) faces the
    # viewer, where a ring spreads its vertices widest.
    marked = [v for v, d in verts.items() if d.get("highlight")]
    turn = -math.pi / 2 - 2 * math.pi * pos[marked[0]] if marked else 0.0
    xy = {}
    for v, d in verts.items():
        th = 2 * math.pi * pos[v] + turn
        i = d["level"]
        xy[v] = (rx[i] * math.cos(th), yc[i] + ry[i] * math.sin(th), math.sin(th))

    counts = defaultdict(int)
    measured = defaultdict(int)
    for d in verts.values():
        counts[d["level"]] += 1
        measured[d["level"]] += (d.get("values") or {}).get(metric["key"]) is not None
    ell = spec["ell"]
    for i in range(depth + 1):
        top, bot = yc[i] + ry[i] + 0.32, yc[i] - ry[i] - 0.18
        ax.add_patch(Rectangle((-width, bot), 2 * width, top - bot, color=BAND, zorder=0, lw=0))
        ax.add_patch(Ellipse((0, yc[i]), 2 * rx[i], 2 * ry[i], fill=False, ec="#dfe3ea", lw=0.8, ls=":", zorder=0.5))
        name = "crater = floor" if depth == 0 else ("crater" if i == 0 else ("floor" if i == depth else ""))
        head = f"Level {i} · {name}" if name else f"Level {i}"
        ax.text(-width + 0.04, top - 0.05, head, fontsize=9, fontweight="bold", color=INK, va="top")
        cond = f"{ell ** i:,}"
        cov = f"{counts[i]} curve{'' if counts[i] == 1 else 's'}" + ("" if measured[i] == counts[i] else f", {measured[i]} measured")
        ax.text(-width + 0.04, top - 0.20, f"{cov} · conductor {cond}", fontsize=7.5, color=MUTED, va="top")

    for (v, w), status in edges.items():
        if v not in xy or w not in xy:
            continue
        (x0, y0, s0), (x1, y1, s1) = xy[v], xy[w]
        ax.plot([x0, x1], [y0, y1], color=EDGE, lw=0.8, ls="--" if status != "verified" else "-",
                zorder=1 + 0.1 * (2 - s0 - s1))

    vals = [(d.get("values") or {}).get(metric["key"]) for d in verts.values()]
    vals = [x for x in vals if x is not None]
    lo, hi = (min(vals), max(vals)) if vals else (0.0, 1.0)
    if hi == lo:
        hi = lo + 1.0
    cmap = plt.get_cmap(CMAP)
    # Label every measured vertex (unless there are very many), and the
    # unmeasured ones too on small walks.
    n_meas = sum(1 for d in verts.values() if (d.get("values") or {}).get(metric["key"]) is not None)
    highlights = []
    for v, (x, y, s) in xy.items():
        d = verts[v]
        val = (d.get("values") or {}).get(metric["key"])
        z = 3 - s  # front of a ring over its back
        if val is not None:
            ax.scatter(x, y, s=95, color=cmap((val - lo) / (hi - lo)), edgecolor="#3a3a3a", lw=0.6, zorder=z)
        else:
            ax.scatter(x, y, s=16, color="#e2e2e2", edgecolor="#9a9a9a", lw=0.5, zorder=z)
        if (val is not None and n_meas <= 120) or len(verts) <= 60:
            ox, oy = x / (rx[d["level"]] or 1), s
            ax.annotate(d.get("label", str(v)), (x, y), textcoords="offset points",
                        xytext=(9 * ox, 9 * oy + (2 if oy >= 0 else -2)), ha="center",
                        va="bottom" if oy >= 0 else "top", fontsize=6.5, color=INK, zorder=6)
        if d.get("highlight"):
            ax.scatter(x, y, s=260, facecolor="none", edgecolor=RING, lw=1.8, zorder=7)
            highlights.append(f"{d.get('label', v)} = {d['highlight']}")

    ax.set_xlim(-width, width)
    ax.set_ylim(yc[-1] - ry[-1] - 0.25, ry[0] + 0.4)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(metric.get("title", metric["key"]), fontsize=12, color=INK, pad=6)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(lo, hi))
    cb = plt.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.045, pad=0.03, shrink=0.6)
    cb.set_ticks([lo, (lo + hi) / 2, hi])
    cb.set_ticklabels([_fmt(lo), _fmt((lo + hi) / 2), _fmt(hi)])
    cb.ax.tick_params(labelsize=8)
    across = "both levels" if depth == 1 else "all levels"
    cb.set_label(f"{metric.get('label', metric['key'])} · one observed-range scale across {across}", fontsize=9, color=INK)
    return highlights


def render(spec, out_stem):
    spec = _v2(spec)
    depth = spec["depth"]
    verts = {v["id"]: v for v in spec["vertices"] if v.get("level") is not None}
    pos, _ = _positions(verts)
    edges = {}
    for v, d in verts.items():
        for w in d.get("nbrs", []):
            if w in verts and v != w:
                edges[(min(v, w), max(v, w))] = "verified"
    for e in spec.get("edges", []):
        edges[(min(e["a"], e["b"]), max(e["a"], e["b"]))] = e.get("status", "verified")

    metrics = spec["metrics"]
    rx, ry, yc = _geometry(depth)
    # Axes are drawn to scale (equal aspect): size the figure to the data,
    # plus room for the title, the colour bar and the footnote.
    span_x, span_y = 2 * (rx[-1] + 0.55), ry[0] + 0.4 - (yc[-1] - ry[-1] - 0.25)
    h = 9.0 * span_y / span_x + 2.4
    fig, axs = plt.subplots(1, len(metrics), figsize=(10 * len(metrics), h), squeeze=False)
    highlights = []
    for ax, m in zip(axs[0], metrics):
        highlights = _panel(ax, spec, m, verts, pos, edges, depth)
    title = spec.get("title", f"{spec['ell']}-isogeny volcano, depth {depth}")
    sub = spec.get("subtitle", "")
    if highlights:
        sub += ("  ·  " if sub else "") + "ringed blue: " + ", ".join(highlights)
    fig.suptitle(title + (f"\n{sub}" if sub else ""), fontsize=15, color=INK)
    foot = spec.get("footnote")
    if not foot:
        per = " above ".join(f"{sum(1 for d in verts.values() if d['level'] == i)}" for i in range(depth + 1))
        foot = (f"{per} curves, crater to floor; grey: not measured; dashed: unverified edge. "
                "Each panel has its own scale.")
    fig.text(0.5, 0.01, foot, ha="center", fontsize=9, color=MUTED)
    fig.tight_layout(rect=(0, 0.03, 1, 1 - 1.0 / h))
    fig.savefig(f"{out_stem}.svg", bbox_inches="tight")
    fig.savefig(f"{out_stem}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    render(json.loads(open(sys.argv[1]).read()), sys.argv[2])
    print("wrote", sys.argv[2] + ".svg", sys.argv[2] + ".png")
