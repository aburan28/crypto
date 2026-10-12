# Audit rendering and operational checks

The SVG was rendered from its source at 150 dpi and the five-page audit PDF
was inspected page by page. The final rendering fits its text, tables and
diagram within the page margins; citations are clickable. Rendering only
lays out existing written proofs and native receipts.

The initial ImageMagick SVG conversion failed because its `rsvg-convert`
delegate was absent. The retained `diagram_render.stderr` records that
failure. The final layout uses CairoSVG installed in a scratch directory;
its source renderer and installation log are retained. A harmless fontconfig
initialization warning is recorded by that renderer.

The first local CI-guard execution used a reduced toolchain PATH and could
not locate the already-installed `rg`. Its `ci_guards.attempt1.log` is
retained. The final execution inherits PATH while prepending the pinned
Lean toolchain, and checks that `rg` exists before the guard runs.

An upstream HTTP 503 tool-transport interruption was reported during the
audit handoff after the proof commands had completed. It does not change
the successful native process exit statuses, Lean build, axiom audit or
fresh kernel replay. No mathematical failure is inferred from that event.

The workflow was parsed and its declared triggers, read-only repository
permission and fresh-checker arguments checked locally. Its source/axiom
guards were exercised locally. Remote CI was not launched or monitored.
