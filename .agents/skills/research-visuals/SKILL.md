---
name: research-visuals
description: Produce and refresh evidence-linked diagrams, graphs, reports, and PDFs for every experiment run (ICMS, ecbench, autolab, benchmarks), every PR that lands a result, and every search for isogenies, curves, scalar rules, endomorphisms, or related ECDLP mechanisms.
---

# Research visuals and reports

Use this skill at three points (AGENTS.md "Research searches must leave
visual reports"):

- **every experiment run**: an ICMS or ecbench session, an autolab or
  tournament round, a benchmark, or any other research run ("Every experiment
  run" below);
- **every PR that lands a result**: a comparison, claim, decision, or
  scoreboard change ("Every result PR" below);
- **every substantive search round**: a literature, theory, code, or
  experiment search for new isogenies, curves, scalar rules, endomorphisms,
  or nearby ECDLP ideas ("Deliverables for each search round" below).

Read `AGENTS.md` first. Keep the repository's native execution, identity,
measurement, PR, and claim gates in force. Python here only lays out what was
recorded (plotting recorded values, rendering, PDF inspection); it never
computes a result (AGENTS.md "Implementation language: no Python").

## Every experiment run

Each run leaves a descriptive run report with a graph and a PDF, in the same
PR as the run's records; failed and timed-out runs included.

- **Where.** Beside the records, never inside a directory a harness seals or
  audits: ICMS sessions are add-only and their `exec/` directories hold
  exactly the files the record hashes, and ecbench CI re-audits committed
  sessions. Use `docs/ic/measurement/reports/<session-directory-name>/` for
  an ICMS session and `research/<study>/reports/run-<YYYYMMDDTHHMMZ>/` (UTC
  start) for an ecbench session, autolab round, benchmark, or other run.
  Write-once: a later run writes its own report.
- **What it says.** The spec or command, its id and commit, the host and the
  isolation level earned, arms, seeds and workloads, runs requested and
  completed with each failure's recorded reason, the measured values with
  their unit, sample size and the spread or interval the harness recorded,
  references and controls, and the paths of the raw records. Label every
  number measured or derived, and show how a derived one was computed. State
  missing data as missing.
- **Observations only.** No speedup, win, or significance language beyond
  what the harness's own comparison admitted, with its scope; wall time only
  where the harness admitted it. Interpretation belongs to the result PR.
- **Graph.** At least one quantitative graph of the recorded values (arms
  side by side, per-round values, a quantity against a parameter), with
  units, sample size, uncertainty where recorded, and the session or run id
  in the caption. Plot recorded values only: no new computation or reruns. A
  run with fewer than two comparable values gets a diagram of what ran
  instead, and says why. Keep the graph's source or data beside its SVG.
- **PDF.** Build it with the graph embedded ("Rendering" below) and inspect it
  before committing.

## Every result PR

A PR that lands a comparison, claim, decision, or scoreboard change carries a
results report in the study directory (the README or REPORT the PR already
needs, for example ecbench's "Land it" README) with the graphs its conclusion
rests on and a PDF: measured values with units, sample sizes and intervals
against the reference, threshold, or incumbent, with session and run ids, and
the scope and limits of the conclusion. Link the run reports it builds on.
Refresh every affected canonical graph and its rendered copies in the same PR
("Refresh existing graphs in the same change" below; AGENTS.md section 7), and
list the graphs checked and left unchanged.

## Rendering

Use an installed renderer: Typst, Pandoc with a PDF engine, or LaTeX. Where
none is installed (a Claude Code cloud container has none), install Typst as a
Python wheel in a scratch virtualenv outside the repository; it is document
rendering, not research code:

```sh
python3 -m venv /tmp/report-venv
/tmp/report-venv/bin/pip install -q typst matplotlib
/tmp/report-venv/bin/python -c "import typst; typst.compile('report.typ', output='report.pdf')"
```

Write the report source in Typst when Typst renders it (Markdown needs
Pandoc). Draw graphs as SVG from the recorded values (matplotlib, Graphviz, a
native renderer, or hand-written SVG) and embed the SVG. Open the PDF and read
every page: it opens, the graph is legible and labeled, and ids and citations
are present. Never commit a PDF you have not inspected or present an older
PDF as current. If no renderer can be installed, commit the complete source
and graphs and state the exact blocker in the report and the PR.

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
   document tool such as Pandoc, Typst, or LaTeX ("Rendering" above). Keep
   source and PDF together.
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
