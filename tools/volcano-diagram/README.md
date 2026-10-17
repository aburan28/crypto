# Standard isogeny-volcano heat map

`volcano_diagram.py` is the repository's one way to draw an ℓ-isogeny
volcano. Every study that walks a volcano, or puts a measurement on its
curves, uses it, so all such figures read the same way (`research-visuals`
skill, "Isogeny volcanoes: the standard diagram").

![example](example/volcano.png)

What it shows:

- **Strata.** Each level is a shaded band, crater on top and floor at the
  bottom. The band is labelled `Level i · crater` / `Level i` /
  `Level i · floor`, with its curve count, how many were measured if not
  all, and its conductor relative to the crater, ℓ^i.
- **Rings.** Inside its band, each level is an elliptical ring seen in
  perspective, each ring wider than the one above, so the figure reads as
  a volcano. Crater vertices go around their ring in cycle order. Each
  lower vertex sits at the angle of the curves it leads down to, under its
  parent, and the front half of a ring is drawn over its back.
- **Panels.** One panel per metric, normally mean operations beside mean
  time. Each panel colours measured vertices on one observed-range YlOrRd
  scale shared by every level, with a colour bar printing min, mid and
  max. Vertices without a value are small and grey, so it's visible what
  was not measured.
- **Edges.** Verified edges are solid grey and unverified ones are dashed.
- **Marks.** A vertex with `highlight` (e.g. the walk's root) is ringed blue
  and named in the subtitle. The title says what was measured; the
  subtitle carries the run or session id; the footnote states the level
  counts or anything else the colours leave out.
- **Output.** `OUT.svg` (vector, for reports and PDFs) and `OUT.png` at 150
  dpi (for `benchmarks/`), from the same figure.

## Input: `volcano-diagram/v2`

```json
{
  "schema": "volcano-diagram/v2",
  "ell": 2,
  "depth": 3,
  "title": "Measured 2-isogeny volcano heat map · p = 146833024135393 · 100 paired walk seeds per curve",
  "subtitle": "run-20261017T0049Z",
  "footnote": "optional; defaults to the per-level counts",
  "metrics": [
    {"key": "ops",  "title": "native rho · mean iterations", "label": "group operations per solve"},
    {"key": "time", "title": "native rho · mean time",       "label": "ms per solve"}
  ],
  "vertices": [
    {"id": 0, "label": "v0", "level": 3, "nbrs": [1], "values": {"ops": 497201.0, "time": 371.5}, "highlight": "root"},
    {"id": 1, "label": "v1", "level": 2, "nbrs": [0, 4], "values": {}}
  ],
  "edges": [{"a": 0, "b": 1, "status": "verified"}]
}
```

- `level`: 0 = crater … `depth` = floor. Vertices with `"level": null` are
  left out.
- `nbrs`: the vertex's neighbours. Edges drawn from it count as verified.
- `values`: one entry per metric `key`; a missing key means not measured.
- `edges` (optional): overrides an edge's status (`verified` or
  `unverified`).
- `label`: a short tag on the drawing. Keep each vertex's full curve id in
  the study's records and map label → id in the report, since a volcano
  has too many vertices to print full ids on the figure.
- `volcano-diagram/v1` input (one `value` per vertex, a `value_label`) is
  still read and drawn as a single panel.

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
"Implementation language: no Python"). The users are `research/isogeny-volcano-rho/plot.py` (from a `volcano-rho`
`results.json`) and `plot_ca.py` (from a `ca volcano` run).
