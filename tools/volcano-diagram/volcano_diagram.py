"""The repository's standard isogeny-volcano diagram (`volcano-diagram/v1`).

Lays out a recorded ℓ-isogeny volcano as horizontal strata (crater on top,
floor at the bottom), hangs every lower vertex under its parent so each
crater vertex's tree reads left to right, and colours the vertices that
carry a measured value on one shared scale. Vertices without a value are
drawn small and grey, so coverage is visible. Plotting only: every number
it draws comes from the input file (AGENTS.md "Implementation language").

usage:
    python volcano_diagram.py INPUT.json OUT_STEM     # writes OUT_STEM.svg and OUT_STEM.png (150 dpi)

or from Python:
    from volcano_diagram import render
    render(spec, "out/volcano")

Input (`volcano-diagram/v1`, see README.md):
    {
      "schema": "volcano-diagram/v1",
      "ell": 3, "depth": 3,
      "title": "3-isogeny volcano, depth 3, p = ...",
      "subtitle": "run-...",                      # run or session id
      "value_label": "mean rho iterations",       # colour-bar label; omit if no values
      "vertices": [
        {"id": 0, "label": "v0", "level": 3, "nbrs": [1],
         "value": 658723.0,                       # optional: measured value
         "highlight": "root"}                     # optional: circled red, named in the title
      ],
      "edges": [{"a": 0, "b": 1, "status": "unverified"}]   # optional; default: nbrs, verified
    }
Levels count down from the crater: 0 = crater, depth = floor. Vertices with
`level: null` are left out. Unverified edges are drawn dashed.
"""

import json
import sys
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SCHEMA = "volcano-diagram/v1"
CMAP = "viridis"
EDGE = "#bbbbbb"


def _layout(verts, depth):
    """x positions: leaves consecutive, parents over their children's mean."""
    crater = sorted(v for v, d in verts.items() if d["level"] == 0)
    children = defaultdict(list)
    for v, d in sorted(verts.items()):
        if d["level"] > 0:
            ups = sorted(w for w in d["nbrs"] if w in verts and verts[w]["level"] == d["level"] - 1)
            if ups:
                children[ups[0]].append(v)
    xpos, cursor = {}, [0.0]

    def place(v):
        if not children[v]:
            xpos[v] = cursor[0]
            cursor[0] += 1.0
        else:
            for c in children[v]:
                place(c)
            xpos[v] = sum(xpos[c] for c in children[v]) / len(children[v])

    for c in crater:
        place(c)
        cursor[0] += 0.6
    return xpos, cursor[0]


def render(spec, out_stem):
    if spec.get("schema") != SCHEMA:
        raise ValueError(f"expected schema {SCHEMA}")
    depth = spec["depth"]
    verts = {v["id"]: v for v in spec["vertices"] if v.get("level") is not None}
    xpos, width = _layout(verts, depth)

    edges = {}
    for v, d in verts.items():
        for w in d.get("nbrs", []):
            if w in verts and v != w:
                edges[(min(v, w), max(v, w))] = "verified"
    for e in spec.get("edges", []):
        edges[(min(e["a"], e["b"]), max(e["a"], e["b"]))] = e.get("status", "verified")

    fig, ax = plt.subplots(figsize=(min(18, max(8, width * 0.12)), 4.8))
    for (v, w), status in edges.items():
        if v in xpos and w in xpos:
            ax.plot([xpos[v], xpos[w]], [-verts[v]["level"], -verts[w]["level"]], color=EDGE, lw=0.7,
                    ls="--" if status != "verified" else "-", zorder=1)

    valued = [d["value"] for d in verts.values() if d.get("value") is not None]
    lo, hi = (min(valued), max(valued)) if valued else (0.0, 1.0)
    cmap = plt.get_cmap(CMAP)
    highlights = []
    for v, x in xpos.items():
        d = verts[v]
        y = -d["level"]
        if d.get("value") is not None:
            ax.scatter(x, y, s=70, color=cmap((d["value"] - lo) / (hi - lo or 1)), edgecolor="k", zorder=3)
            floor = d["level"] == depth
            ax.annotate(d.get("label", str(v)), (x, y), textcoords="offset points", xytext=(0, -9 if floor else 7),
                        ha="center", va="top" if floor else "bottom", fontsize=6, rotation=90 if floor else 0)
        else:
            ax.scatter(x, y, s=12, color="#dddddd", edgecolor="#999999", zorder=2)
        if d.get("highlight"):
            ax.scatter(x, y, s=200, facecolor="none", edgecolor="red", lw=1.5, zorder=4)
            highlights.append(f"{d.get('label', v)} = {d['highlight']}")

    if valued:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(lo, hi))
        fig.colorbar(sm, ax=ax, label=spec.get("value_label", "value"))
    names = ["crater (0)"] + [f"level {lv}" for lv in range(1, depth)] + ([f"floor ({depth})"] if depth else [])
    ax.set_yticks([-lv for lv in range(depth + 1)], names)
    ax.set_xticks([])
    title = spec.get("title", f"{spec['ell']}-isogeny volcano, depth {depth}")
    if highlights:
        title += f"  (circled red: {', '.join(highlights)})"
    if spec.get("subtitle"):
        title += f"\n{spec['subtitle']}"
    ax.set_title(title, fontsize=9)
    fig.savefig(f"{out_stem}.svg", bbox_inches="tight")
    fig.savefig(f"{out_stem}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    render(json.loads(open(sys.argv[1]).read()), sys.argv[2])
    print("wrote", sys.argv[2] + ".svg", sys.argv[2] + ".png")
