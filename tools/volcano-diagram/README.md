# Standard isogeny-volcano diagram

`volcano_diagram.py` is the repository's one way to draw an ℓ-isogeny
volcano. Every study that walks a volcano, or puts a measurement on its
curves, uses it, so all such figures read the same way (`research-visuals`
skill, "Isogeny volcanoes: the standard diagram").

![example](example/volcano.svg)

What it shows:

- **Strata.** Horizontal rows by level: crater (0) on top, the floor at the
  bottom, and every level between labelled.
- **Crater.** Crater vertices are laid out in cycle order. An edge between
  non-adjacent crater positions (the one closing a full cycle) arcs above.
- **Trees.** Each lower vertex hangs under its first neighbour one level up.
  Leaves are spaced evenly and parents sit over their children, so each
  crater vertex's descending tree reads left to right.
- **Values.** Vertices with a measured `value` are drawn large and coloured
  on one shared viridis scale with a labelled colour bar. Vertices without
  one are drawn small and grey, so it's visible what was not measured.
- **Edges.** Verified edges are solid grey and unverified ones are dashed.
- **Marks.** A vertex with `highlight` (e.g. the walk's root) is circled
  red and named in the title. The title carries ℓ, the depth and the field;
  the subtitle carries the run or session id.
- **Output.** `OUT.svg` (vector, for reports and PDFs) and `OUT.png` at 150
  dpi (for `benchmarks/`), from the same figure.

## Input: `volcano-diagram/v1`

```json
{
  "schema": "volcano-diagram/v1",
  "ell": 3,
  "depth": 3,
  "title": "3-isogeny volcano, depth 3, p = 1090993606999",
  "subtitle": "run-20261016T2302Z",
  "value_label": "mean rho iterations (measured curves)",
  "vertices": [
    {"id": 0, "label": "v0", "level": 3, "nbrs": [1], "value": 658723.0, "highlight": "root"},
    {"id": 1, "label": "v1", "level": 2, "nbrs": [0, 4], "value": null}
  ],
  "edges": [{"a": 0, "b": 1, "status": "verified"}]
}
```

- `level`: 0 = crater … `depth` = floor. Vertices with `"level": null` are
  left out.
- `nbrs`: the vertex's neighbours. Edges drawn from it count as verified.
- `edges` (optional): overrides an edge's status (`verified` or
  `unverified`).
- `label`: a short tag on the drawing. Keep each vertex's full curve id in
  the study's records and map label → id in the report, since a volcano
  has too many vertices to print full ids on the figure.

Save the input next to the figure (`<name>.diagram.json`), so the figure can
be re-rendered from the input alone.

## Use

```sh
python tools/volcano-diagram/volcano_diagram.py run/volcano.diagram.json run/volcano
```

or from a study's plotting script:

```python
sys.path.insert(0, "<repo>/tools/volcano-diagram")
from volcano_diagram import render
render(spec, "run/volcano")
```

It needs matplotlib only, and plots recorded values only (AGENTS.md
"Implementation language: no Python"). The first user is
`research/isogeny-volcano-rho/plot.py`, which builds the input from a
`volcano-rho` `results.json`.
