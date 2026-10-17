// Results report for research/generic_dlp_memory_parallelism_20261017
// (research-visuals skill, "Every result PR").  Every number is read from
// the run's summary.json, written by examples/generic_dlp_analyze.rs.
#let run = "reports/run-20261017T0201Z/"
#let s = json(run + "summary.json")
#set page(paper: "a4", margin: 1.7cm, numbering: "1")
#set text(size: 9.5pt)
#set heading(numbering: "1.")
#let f3(x) = if x == none { "–" } else { str(calc.round(x, digits: 3)) }
#let f4(x) = if x == none { "–" } else { str(calc.round(x, digits: 4)) }
#let cfg(part, bits, kernel) = s.configs.find(c => c.part == part and c.bits == bits and c.kernel == kernel)

#align(center)[
  #text(15pt, weight: "bold")[Memory-level parallelism in BSGS, inversion-free rho start-up]

  Study `research/generic_dlp_memory_parallelism_20261017` · run `run-20261017T0201Z` · 2026-10-17 · PR \#2048
]

*Class: engineering* (AGENTS.md §3).
- Two changes to the single-word generic-DLP engines lower one-target wall
  time.
- Every counted output is identical between baseline and candidate: answers,
  steps, table entries, candidates, distinguished points, and the field
  operations per visited point.
- `S` and every IC/rho ratio are unchanged.
- Host: one x86-64 Xeon (Sapphire/Emerald Rapids class) at 2.1 GHz, a 4-vCPU
  KVM guest, isolation L0 with contended runs left out.

= What changed

+ *`bsgs_fast`: table slots are prefetched a round before they are probed.*
  - The baby-step table outgrows L2 past about `2^36`.
  - At load factor ½ the probe's "empty?" branch depends on the missing
    load, so successive misses serialise.
  - The kernels now hint each visited point's `x` to the table before the
    visit. For `paired` and `co-Z`, a round is computed and hinted in one
    pass and visited in a second, in the original order.
  - Visits, early exit, charged steps and field operations are unchanged;
    a test pins all four.
+ *`prime_fast`: Jacobian scalar multiplication, Straus–Shamir `combine`.*
  - `rho_parallel` builds its 1024-entry `r`-adding table and every walk
    start with `combine(g, q, u, v)`. That was two affine double-and-add
    ladders paying a binary-GCD inversion per group operation.
  - None of this work is in `total_steps`, and all of it is target-dependent.
  - It now costs two inversions per call. The results equal the affine
    ladder's on every point of four small curves and on the ladder curves.

*Tried and rejected.* `MADV_HUGEPAGE` on the table.
- In exploratory runs it sped up the search phase but made the build phase
  erratic: THP defrag is `madvise` here, so the first touch compacts
  synchronously.
- It was not shipped.

= Method

- *Arms.* A is `main` at `b1bf5b893`; A′ is the same binary, as the A/A
  control; B is commit `c699775dd`.
- *Rounds.* 7 rotated rounds, 5 at 52 bits.
- *Workload.* One target per process, a cold full-group solve.
- *Isolation.* Every run is under `tools/isolated_bench.py`.
- *Rule.* "#s.rule".
- *Contention.* #s.configs.map(c => c.contended.A + c.contended.A2 + c.contended.B).sum() of 438 timed runs were contended and are left out, as are #s.pari_contended_left_out of 30 PARI runs.
- *Pinned outputs.* Identical in every configuration: #if s.configs.all(c => c.pinned_identical) [yes] else [*no*].
- *Field-operation rows.* #s.opcount_rows_compared rows compared, #if s.opcount_identical [all identical] else [*mismatch*].
- *Curves.* They are in `inputs/curves.json`, with `#E` checked by PARI/GP:
  - BSGS uses `y² = x³ + 7` (W11): `icv1-fp40-t1502959-42c6aada`,
    `icv1-fp44-tm1269095-f56aca35`, `icv1-fp48-t14984575-af5d5615`,
    `icv1-fp52-tm131765273-051185e0`;
  - rho uses `a = −3`: `icv1-fp32-t7319-90da43e9`, `icv1-fp36-t454193-44a20818`,
    `icv1-fp40-tm924137-0ecaaeb9`, `icv1-fp44-tm6431881-43acc95d`,
    `icv1-fp48-tm15186475-c683ebfc`, `icv1-fp52-t113543683-af85fb1e`.

= Results

== BSGS, one target, full group, negation map

Times are in seconds; ratios are A/B. #text(fill: gray)[The per-round ranges are min–max.]

#table(
  columns: 8, stroke: 0.3pt, inset: 3pt,
  [threads], [bits], [kernel], [A `total_s`], [B `total_s`], [A/B], [A/B per round], [verdict],
  ..s.configs.filter(c => c.part.starts-with("bsgs")).map(c => {
    let t = c.timers.total_s
    (str(c.threads), str(c.bits), c.kernel, f4(t.A.median), f4(t.B.median), f3(t.ratio_of_medians_A_over_B),
     "[" + f3(t.per_round_A_over_B.min) + ", " + f3(t.per_round_A_over_B.max) + "]", t.verdict)
  }).flatten()
)

The giant (search) phase moved by about as much as the build. The frozen
targets have small logarithms, so `total_s` is mostly the table build; giant
steps are pinned identical, so a `search_s` ratio is a per-step ratio. The
search ratios, from the run report:

- paired: #f3(cfg("bsgs1", 44, "paired").timers.search_s.ratio_of_medians_A_over_B) at 44 bits,
  #f3(cfg("bsgs1", 48, "paired").timers.search_s.ratio_of_medians_A_over_B) at 48,
  #f3(cfg("bsgs1", 52, "paired").timers.search_s.ratio_of_medians_A_over_B) at 52;
- chain: #f3(cfg("bsgs1", 48, "chain").timers.search_s.ratio_of_medians_A_over_B) at 48 bits.

#figure(image(run + "fig_bsgs.svg", width: 100%), caption: [BSGS medians and A/B with per-round ranges against the A/A′ noise band (run-20261017T0201Z).])

== Rho (`prime_fast::rho_parallel`), one target, negation map

#table(
  columns: 8, stroke: 0.3pt, inset: 3pt,
  [threads], [bits], [A `solve_s`], [B `solve_s`], [A/B], [A/B per round], [verdict], [start-up alone, A/B],
  ..s.configs.filter(c => c.part.starts-with("rho")).map(c => {
    let t = c.timers.solve_s
    (str(c.threads), str(c.bits), f4(t.A.median), f4(t.B.median), f3(t.ratio_of_medians_A_over_B),
     "[" + f3(t.per_round_A_over_B.min) + ", " + f3(t.per_round_A_over_B.max) + "]", t.verdict,
     f3(c.timers.setup_s.ratio_of_medians_A_over_B))
  }).flatten()
)

#figure(image(run + "fig_rho.svg", width: 100%), caption: [Rho solve and start-up medians, and A/B at 1 and 3 threads (run-20261017T0201Z).])

== External reference: PARI/GP `elllog` (AGENTS.md §15)

PARI/GP 2.15.4, single-threaded, on the same curves and the same targets,
timed by `getabstime()` around `elllog`. Every answer was verified. Which
generic method PARI ran is not observable: its help text and debug domains
print none.

#table(
  columns: 7, stroke: 0.3pt, inset: 3pt,
  [curves], [bits], [PARI median (s)], [A median (s)], [B median (s)], [B / PARI], [A / PARI],
  ..s.external.map(e => (e.role, str(e.bits), f3(e.pari_elllog_s.median), f4(e.native_A_median_s),
    f4(e.native_B_median_s), f4(e.B_over_pari), f4(e.A_over_pari))).flatten()
)

Growth is `log₂` of the median-time ratio per bit; `√n` scaling is 0.5:

#for g in s.external_growth [
  - #g.role curves, #g.from_bits→#g.to_bits bits: PARI #f3(g.pari_log2_per_bit), A #f3(g.native_A_log2_per_bit), B #f3(g.native_B_log2_per_bit).
]

#figure(image(run + "fig_external.svg", width: 100%), caption: [PARI/GP `elllog` against both native arms on the same curves and targets.])

*Scope of the external comparison.*
- At every size measured, both native arms are faster than PARI/GP on this
  host. Over this range the native rho grows more slowly than `√n`, because
  its fixed start-up is amortised, so the 52-bit gap is not a constant to
  extrapolate.
- The fastest published ECDLP implementations (record-scale
  distinguished-point rho on GPUs and FPGAs, at 100 bits and more) were not
  reproduced.
- Sage is not installed, and Magma needs a licence that this host lacks.
- No "fastest" claim is made.

= Interpretation (labelled)

- *BSGS was bound by memory latency, not arithmetic.* Taking the probe's
  miss off the dependency chain saved 23–31% (chain) and 43–50% (paired) of
  one-target time at 44–52 bits.
  - Paired gained more: its arithmetic is 3.8–4.5× smaller (3.84–3.94
    against 14.6–17.9 M-equivalents per point, from this run's op counts),
    so it was the more memory-bound kernel.
  - At three threads, which share one table and one memory system, the gain
    is smaller.
- *Rho at small sizes was dominated by its uncounted start-up.*
  - That start-up is now 6.8–16.1× cheaper.
  - The solve gain shrinks as the `√n` walk grows. At 44 and 52 bits (one
    thread) it is within noise.

= Limits

- No change to `S`, to a counted cost, or to any IC/rho comparison. No
  method ranking; the §5a band does not apply.
- One host class only. Arm64, GPUs, other x86 parts and bare metal were not
  measured.
- One frozen target per curve, as the protocol fixes; the targets have small
  logarithms. The 40-bit BSGS rows are no measurable change by the rule.

= Audit and next steps

The audit behind this study covers BSGS, rho, kangaroo, Gaudry–Schost,
Pohlig–Hellman and the `ecc/` stacks. It found 26 findings in all, with how
each is known and its status, in `AUDIT_FINDINGS.md`.

- *Fixed in a separate PR (pending):*
  - the `pollard_collab` negation-map 2-cycle escape, which never fired;
  - Pohlig–Hellman's `O(q)` digits, unbounded factorisation and unchecked
    answer.
- *Open:*
  - a coarse scalar bit-length timing leak in `Point::scalar_mul_ct`'s
    generic fallback;
  - Gaudry–Schost and collab branch bits overlapping the DP bits;
  - grumpy giants not following Bernstein–Lange;
  - `aut_folded_rho` not folding;
  - several reference-arm cost levers, each needing its own calibration.
- *Next engineering steps, named in the README:*
  - interleaved negation-map BSGS for one target (`≈ 0.943√n` against
    `1.00√n`, derived);
  - a barrier-free single-worker baby phase;
  - coefficient-free rho steps;
  - a Jacobian `bsgs_fast` scalar multiplication.
