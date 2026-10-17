// Run report for run-20261017T0201Z (research-visuals skill, "Every
// experiment run").  Descriptive only: every number is read from
// summary.json (written by examples/generic_dlp_analyze.rs) or
// manifest.json, never computed here.
#let s = json("summary.json")
#let m = json("manifest.json")
#set page(paper: "a4", margin: 1.6cm, numbering: "1")
#set text(size: 9pt)
#set heading(numbering: "1.")
#let f3(x) = if x == none { "–" } else { str(calc.round(x, digits: 3)) }
#let f4(x) = if x == none { "–" } else { str(calc.round(x, digits: 4)) }

= Run report: #m.run

*Study* `research/generic_dlp_memory_parallelism_20261017` ·
*protocol* `PROTOCOL.md` (frozen before timing; amendment 1) ·
*Conductor task* T-1 · *PR* #2048.

This report describes what ran and what was recorded. The study's
`README.md` and `REPORT.pdf` carry the interpretation.

== What ran

- *Baseline A and A′* (the A/A control): `main` at #raw(m.baseline_commit).
- *Candidate B*: #raw(m.candidate_commit) on #raw(m.candidate_branch).
- *Build*: #m.build, #m.rustc.
- *Host*: #m.host.cpu, #m.host.logical_cpus logical CPUs, #m.host.virtualization.
  - Features: #m.host.flags.
  - Caches: #m.host.caches.
  - Kernel #m.host.kernel; transparent huge pages #raw(m.host.transparent_hugepage).
- *Isolation* (as written at run start): #m.isolation.
  - *Amendment*: #m.amendments.join(" ")
  - Isolation level earned: shared KVM guest, so L0 by ecbench's scale.
  - Contended runs are flagged per record and left out of every statistic.
- *External tool*: #m.external_tool.
- *Commands*: #for c in m.commands [#raw(c) \ ]
- *Runs*:
  - BSGS: 1 thread at 40, 44, 48 and 52 bits; 3 threads at 48 bits. Kernels
    `chain` and `paired` (`K = 16`).
  - Rho: 1 and 3 threads at 32–52 bits.
  - 7 rounds (5 at 52 bits), arms rotated `A B A′ / B A′ A / A′ A B`.
  - One target per process.
  - PARI/GP: 3 repetitions per curve.
- *Interruptions* (`interruptions.txt`), none of which reran or overwrote a
  recorded sample:
  1. A lock conflict stopped the script, which then resumed.
  2. The live output was moved out of the repository between samples.
  3. The multi-thread parts started under amendment 1: 3 threads on CPUs
     1–3, because the isolation tool refuses to reserve every CPU.

== Recorded values

Times are in seconds, from each binary's own timer. *A/B* is the ratio of
the medians of A and B. *Per-round* gives the range of the per-round paired
ratios. The *verdict* applies the protocol's rule: "#s.rule".

#let rows = s.configs
#table(
  columns: 11,
  stroke: 0.3pt,
  inset: 2.5pt,
  [part], [bits], [kernel], [timer], [A med.], [B med.], [A/B], [A/B per round], [A/A′ per round], [verdict], [pinned],
  ..rows.map(c => c.timers.pairs().map(((t, v)) => (
    c.part, str(c.bits), c.kernel, t,
    f4(v.A.median), f4(v.B.median), f3(v.ratio_of_medians_A_over_B),
    "[" + f3(v.per_round_A_over_B.min) + ", " + f3(v.per_round_A_over_B.max) + "]",
    "[" + f3(v.per_round_A_over_A2.min) + ", " + f3(v.per_round_A_over_A2.max) + "]",
    v.verdict, if c.pinned_identical { "identical" } else { "MISMATCH" },
  ))).flatten().flatten()
)

*Samples and contention.* Each row below gives the samples per configuration
and the contended runs left out per arm (A, A′, B).

#table(
  columns: 6, stroke: 0.3pt, inset: 2.5pt,
  [part], [bits], [kernel], [samples], [contended (A, A′, B)], [pinned fields],
  ..rows.map(c => (c.part, str(c.bits), c.kernel, str(c.samples),
    str(c.contended.A) + ", " + str(c.contended.A2) + ", " + str(c.contended.B),
    c.pinned_fields.join(", "))).flatten()
)

*Field operations per visited point.* `lit2026_w11_bsgs_opcount` produced
#s.opcount_rows_compared B rows, each compared with its A row. They are
#if s.opcount_identical [identical] else [*not identical*] in
multiplications, squarings, inversions, additions, points visited and slow
rounds.

#figure(image("fig_bsgs.svg", width: 100%), caption: [BSGS, #m.run: one-target `total_s` medians per arm (left) and A/B with the per-round range against the A/A′ range (right). 7 rounds, 5 at 52 bits.])

#figure(image("fig_rho.svg", width: 100%), caption: [Rho, #m.run: `solve_s` and the start-up `combine` calls timed alone (left), and A/B with the per-round range at 1 and 3 threads (right).])

== External reference (PARI/GP `elllog`)

Same curves and targets, single-threaded. Times are medians, in seconds.
PARI is timed by `getabstime()` around `elllog` alone.

#table(
  columns: 8, stroke: 0.3pt, inset: 2.5pt,
  [curve role], [bits], [PARI med. (n)], [A med.], [B med.], [B / PARI], [A / PARI], [PARI verified],
  ..s.external.map(e => (e.role, str(e.bits),
    f3(e.pari_elllog_s.median) + " (" + str(e.pari_elllog_s.values.len()) + ")",
    f4(e.native_A_median_s), f4(e.native_B_median_s), f4(e.B_over_pari), f4(e.A_over_pari),
    if e.pari_verified { "yes" } else { "NO" })).flatten()
)

The native arm is the paired kernel's `total_s` for the BSGS curves and
`solve_s` for the rho curves, both at one thread.

#figure(image("fig_external.svg", width: 100%), caption: [PARI/GP 2.15.4 `elllog` against both native arms on the same curves and targets, #m.run.])

== Raw records

`raw/bsgs1.jsonl`, `raw/bsgs3.jsonl`, `raw/rho1.jsonl`, `raw/rho3.jsonl`,
`raw/pari.jsonl`, `raw/opcount.jsonl` and `raw/isolation.jsonl` (one
isolation record per timed process, joined to its sample by label);
`summary.json`; `manifest.json`; `interruptions.txt`.
