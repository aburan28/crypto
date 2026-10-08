---
name: research-visuals
description: Produce and refresh evidence-linked diagrams, graphs, reports, and PDFs when searching for isogenies, curves, scalar rules, endomorphisms, or related ECDLP mechanisms.
---

# Research visuals and reports

Use this skill for a substantive literature, theory, code, or experiment search
for new isogenies, curves, scalar rules, endomorphisms, or nearby ECDLP ideas.
Read `AGENTS.md` first. Keep the repository's native execution, identity,
measurement, PR, and claim gates in force.

## Deliverables for each search round

1. Write a dated Markdown report in the relevant `research/` or `docs/`
   study directory. State the question, scope, method, exact curve and subgroup
   identities, sources and run IDs, findings (including negatives), unresolved
   checks, and whether each statement is proposed, derived, independently
   verified, or measured. For a scalar rule, give its domain, preconditions,
   formula, and verification or counterexample; for an isogeny, give source and
   target IDs, degree, map or kernel evidence, and subgroup transport status.
2. Add at least one labeled diagram that explains the search result, plus a
   quantitative graph when the result has data to compare.
   Keep editable source beside its rendered SVG or other vector output. Put
   full curve IDs on isogeny vertices, degree and direction on their edges,
   and explicit status on unverified links. Include units, population, and uncertainty for
   quantitative plots. Make the report link the visual and its evidence.
3. Build a PDF of that report with the visual included, using an available
   document tool such as Pandoc, Typst, or LaTeX. Keep source and PDF together.
   Check that the PDF opens and the labels, citations, and graphs are legible.
   If the renderer is unavailable, save complete source and report the build
   blocker; do not present an old PDF as current.
4. If the round claims an algebraic identity (a bilinear identity, a
   summation-polynomial relation, a counting argument another result rests
   on), ship an `identity.certificate/v1` record beside the claim and cite
   its `IDC1h…` id, per `docs/identity-certificates/README.md`. A claim
   without one is prose.

## Refresh existing graphs in the same change

- Search for every affected canonical graph, including
  `docs/index-calculus-scoreboard.html`, `docs/ic/progress-timeline.json`,
  relevant `docs/curves/` diagrams, `docs/performance-gains/`, `figures/`, and
  study-local figures. Follow `AGENTS.md` section 7 for scoreboard data and
  generated copies. Update source data first, then regenerate all rendered
  forms; never hand-edit a generated copy alone.
- Cite frozen evidence for every new value or edge. Preserve prior values as
  historical entries where the graph contract requires them. Mark proposals,
  extrapolations, and unverified routes clearly; do not draw them as proven
  connections or measured speedups. Do not merge distinct curve models by bit
  count or alias.
- If a search finds no new graphable fact, say so in the report and identify
  the graphs checked. A negative result still gets its report, diagram, and
  PDF; a canonical graph changes only if its conclusion or coverage changes.
- Review the report, vector visual, PDF, graph data, and any published copy
  together in the same PR as the finding. Apply existing validation and merge
  gates; a publication alone is not the source of truth.
