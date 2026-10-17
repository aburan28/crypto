---
name: isogeny-volcano-diagrams
description: Render evidence-linked isogeny volcanoes at arbitrary depth, with exact vertex identities, directed maps, and shared-scale operation or time heatmaps.
---

# Isogeny volcano diagrams

Use this skill for a volcano diagram or heatmap, including multiple levels,
multiple surface vertices, horizontal cycles, disconnected components, or
parallel maps. Read the repository AGENTS.md and research-visuals skill first.

This skill renders recorded graph and measurement data. It does not discover
isogenies, prove volcano structure, solve logarithms, or infer missing levels.
Do not generate synthetic measurements to fill a real graph.

## Inputs and identities

Prepare a JSON object with `title`, `run_id`, `caption`, `nodes`, `edges`,
and optional `metric` and `range`. Each node has a unique string `id`, a
nonnegative integer `level`, and its full exact `curve_id`. Optional `details`
retains field, model, subgroup order, conductor, and evidence paths. Each
directed edge has `source`, `target`, positive integer `degree`, and
`verified` boolean; optional `id` and `evidence` retain map identifiers.
Every directed map is an entry, including the dual, loops and parallel maps.

For a heatmap, `metric` specifies `label`, `unit`, `boundary`, and `statistic`.
Each node then has `value` (a finite number or null) and `n` (sample count).
Keep zero distinct from missing. Put min/mean/max, failures and coverage in
`details`; label whether samples are paired targets or independent repetitions.
For comparable panels supply the same explicit `[minimum, maximum]` range.
Never mix cold setup+query with reused-setup online cost in the same scale.
Store large field elements and point orders as decimal strings.

## Render

Run Node on the bundled script, without installing packages:

```
node .agents/skills/isogeny-volcano-diagrams/render.mjs input.json new-output-dir
```

The directory is write-once. Outputs are `volcano.svg`, offline interactive
`index.html`, exact `graph_data.json`, and `render_manifest.json`. The renderer
uses all recorded levels rather than a fixed depth. It retains all vertices
and maps, renders unverified maps as dashed, and displays full curve identities
in tooltips and the interactive inspector. Levels and categorical topology
colors are distinct from a metric heatmap. Hatching means missing, not zero.
A constant metric uses one middle color. Inputs outside a shared explicit
range are rejected instead of silently clipped.

For PDF and a 150 dpi PNG from the same SVG:

```
python3 .agents/skills/isogeny-volcano-diagrams/export.py new-output-dir
```

Python here performs document rendering and inspection only. Requires
Playwright, Chromium and PyMuPDF. Inspect every generated PDF page. For very
wide graphs, retain the vector SVG and make zoomed panels; never remove nodes
or make unreadable labels appear complete. Run `node .../test.mjs` after a
renderer change; the tests exercise depth, surface cycles, loops, parallel
maps, zero/missing values and constant/shared color ranges.

## Evidence and reporting

Keep the data, source, caption, hash manifest and visual beside the research
report. Include the exact run/fixture identities, degrees, direction,
measurement boundary, units, coverage, failures and isolation level. Link
raw evidence without publishing restricted scalar witnesses or RNG tapes.
Do not treat short curve aliases as registry identities or match aliases
between unrelated runs. Do not interpret rho spread as intrinsic hardness
without paired repeats and appropriate controls. A graph of a selected
ell-component is not the whole isogeny class over all degrees.

Follow research-visuals for source-linked reports and inspected PDFs. Timing
results also require the report/graphs/150 dpi PNG copies and index row under
benchmarks/. Keep earlier snapshots immutable; generate a new output directory
when measurements or coverage change. A plot does not promote an experiment
to a verified speedup or complete an unresolved IC recovery requirement.
