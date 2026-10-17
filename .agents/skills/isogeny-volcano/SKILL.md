---
name: isogeny-volcano
description: Walk the ℓ-isogeny volcano (ℓ = 2 or 3) of a prime-field curve — a given one or one searched for by depth and crater shape — carry one ECDLP instance by verified Vélu isogenies to every stratum, measure Pollard rho on each curve (min / mean / max per curve and per level, with isomorphic controls), and land the standard volcano diagram, report, PDF and benchmark entry. Use whenever someone asks to walk an isogeny graph or volcano, compare ECDLP cost across isogenous curves or volcano levels, or run "the volcano thing" on a curve.
---

# Isogeny volcano: walk it, measure ECDLP on every stratum

Read `AGENTS.md` first: the success criterion, §5a, "Research searches
must leave visual reports" and "Performance results also land in
`benchmarks/`". Everything here is native code with fixed seeds. Python only
draws recorded values.

## The tools

| | crypto (this repo) | aburan28/cryptanalysis |
|:--|:--|:--|
| binary | `target/release/volcano-rho` (`src/bin/volcano_rho.rs`) | `build/ca volcano` (`ca_volcano.h`, `src/volcano.c`) |
| curve | `--curve P,A,B [--order N]`, or a search: `--bits --depth --crater any\|split\|inert\|ramified --max-cofactor` | `--p P --a A --b B [--order N]` |
| rho | built-in single-walker 32-adding DP walk, no negation map; or `--engine ca --ca PATH [--ca-arg ...]` | `ca_rho_solve` in-process (negation map on; `--no-negation`, `--threads`, `--dp-bits`, `--r`) |
| output | `--out DIR`: `results.json`, `trials.csv`, `instances.jsonl` | JSON on stdout, table on stderr, `--trials-out CSV`, `--diagram FILE` |
| graphs | `research/isogeny-volcano-rho/plot.py DIR` (writes `volcano.diagram.json` and the SVG/PNGs) | `--diagram` writes `volcano-diagram/v1` input (`plot_ca.py` turns it into v2) |

The two walks were written independently and agree: both find curve A of
`research/isogeny-volcano-rho/` as 80 curves, 2 / 6 / 18 / 54 from crater to
floor. Use one to check the other when a result matters.

Common options: `--ell 2|3`, `--seed S` (curve search, generator, secret
and every walk seed), `--trials T`, `--per-level L` (curves measured per
level, root always included), `--controls C` (isomorphic re-scalings of
the root: same curve, new coordinates).

## 1. Frame it before running

Write down, in the note or PR:

- **Null hypothesis.** Isogenous curves over `F_p` have the same `#E`, and
  an ℓ-isogeny with ℓ ∤ n is an isomorphism on `⟨P⟩`. So generic rho
  should cost the same on every stratum. The controls measure the spread
  a mere change of coordinates gives. A level effect must beat that
  spread.
- **Reference.** Expected group operations are about `√(πn/2)` without a
  negation map and `√(πn/4)` with one. An r-adding walk adds roughly
  `1/√(1 − 1/r)` (×1.016 at r = 32).
- **Scope.** This is rho against rho. It is not an index-calculus
  improvement, has no §5a band, and a relabelled or re-scaled curve is not
  a speedup (CLAUDE.md success criterion).

## 2. Choose the curve

- **Given curve:** `volcano-rho --curve P,A,B` or `ca volcano --p P --a A
  --b B`. Requires a prime `3 < P < 2^62`. Points are counted unless
  `--order` is given; a wrong order is refused.
- **Depth:** it is `v_ℓ(f)` with `t² − 4p = f²·D_K`, printed by both
  tools. Depth 0 means the walk covers crater curves only: try the other
  ℓ, or search for another curve.
- **Search:** `volcano-rho --bits 40 --depth 3 [--crater split]`. A split
  crater is a cycle and gives many crater curves; ramified gives two;
  inert gives one.
- **ℓ ∤ n:** ℓ must not divide the subgroup order n; both tools refuse it.
- **Cost:** one solve costs about √n group operations. n ≈ 2^38 is under
  a second per solve. Size `--trials × curves` to the time you have.

## 3. Walk and verify

- **Edge checks:** every edge is a Vélu isogeny from a rational root of
  the ℓ-division polynomial. Each one is checked: the image lies on the
  codomain, `φ(P+Q) = φ(P) + φ(Q)`, and `n·φ(P) = O`. Report
  `edges_verified`, which must equal the number of non-root curves walked.
- **Levels:** a curve's level is `depth − (distance to the nearest curve
  with one rational ℓ-kernel)`. Both tools refuse a walk where an edge
  skips a level.
- **Special j:** codomains with j = 0 or 1728 are skipped and counted.
- **Supersingular:** a supersingular curve (t = 0) is warned about,
  because its levels need not follow an ordinary volcano.
- **Truncation:** a split crater can be long, so the walk is truncated
  at its radius (default `2·depth + 2`). Say "within the walk radius"
  when it is.

## 4. Measure

- **Arms:** use the root, `--per-level 10`, `--controls 3`, and
  `--trials 100` for a result (`--trials 5` for a smoke run).
- **Seeds:** trial `s` uses the same seed on every curve. Both tools
  interleave curves within each trial so drift spreads evenly.
- **Correctness:** every solve must return the planted `k`. Both tools
  exit non-zero otherwise.
- **Small groups:** on tiny groups `ca` reports a fixed setup cost (about
  400 operations at n ≈ 4,000). The built-in walk lowers its DP rate there
  (`min(8, log₂n/2 − 6)` bits).
- **Timing:** wall time on a shared host is indicative only; cite group
  operations. For a timed claim use `ecbench-measure`/`isolated_bench`.

## 5. Read the results honestly

- **Min and max:** a single solve ranges over ~250×, so min and max are
  the tails of a wide distribution. Compare means with their 95% CIs,
  and level against level and against the control.
- **Many intervals:** with several levels you are looking at several
  intervals. One CI just missing the reference is not a finding.
- **The instance:** each run fixes one (P, k) per curve, so its spread is
  the walk's own. `volcano-rho --instance-study` separates the instance from
  the walk: it redraws k, the generator, both, or a small k. On curve C,
  5,000 solves found no difference (`run-20261017T0415Z-instance`). Cite
  that run rather than re-asking whether P or k matters.
- **Replay:** re-run the same command with `--trials 5` and compare
  iteration counts with the recorded `trials.csv`; they must be
  identical. Cite the match (e.g. 160/160).

## 6. Visuals and report (research-visuals)

- **Run directory:** put each run in
  `research/<study>/runs/run-<YYYYMMDDTHHMMZ>/`. Write-once: a new run gets
  a new directory.
- **Volcano heat map:** draw it with the standard renderer,
  `tools/volcano-diagram/volcano_diagram.py`, from `volcano.diagram.json`
  (`volcano-diagram/v2`): banded strata, rings in perspective, and a mean
  operations panel beside a mean time panel.
  Measured curves are coloured by mean operations; keep the input beside
  the SVG and PNG.
- **Graphs:** also add per-curve (min–max bar, mean ± CI, every trial) and
  per-level graphs (`plot.py` makes all three).
- **Report:** write a README covering:
  - the question;
  - p, a, b, #E, n, t, D, depth, crater shape and k;
  - the walk's size and per-level counts;
  - tables of min / mean ± CI / median / max per level and per curve;
  - the reference value;
  - the reading, labelled measured or derived;
  - the limits;
  - the reproduce command.
- **PDF:** render with Typst (pandoc → typst, then `typst.compile`), and
  read every page. A line starting with `/` breaks Typst, so reword it.
- **Benchmarks:** copy the PDF, SVGs and 150-dpi PNGs to
  `benchmarks/<YYYYMMDD>-<slug>/` and add a row to `benchmarks/README.md`
  in the same PR.

## 7. Land it

- **PR:** open it ready for review and stay with it (CLAUDE.md "Default
  workflow").
- **Size:** the evidence files (results, trials, SVGs) make the diff
  large. That is expected, and review bots may decline it.
- **Code changes:** a change to either binary must replay the recorded
  runs before it lands.
