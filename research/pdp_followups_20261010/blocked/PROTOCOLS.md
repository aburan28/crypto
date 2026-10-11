# Items 1, 2, 4, 7: preregistered protocols, blocked in this environment

Frozen 2026-10-11.  Each item below has its inputs, falsifier and decision
rule fixed now so that a session with the missing resource can run it
without redesign.  The blockers are environmental, not scientific: the
network policy of this cloud environment answers `403` to every
`github.com` repository not attached to the session (checked for
`mtrimoska/WDSat`, `meelgroup/bosphorus`, `meelgroup/linpb` and the
Xorricane repository on 2026-10-10), `crates.io` likewise, and no GPU is
attached (`nvidia-smi` absent).  PyPI and the Ubuntu mirror are reachable,
which is how cvc5, CryptoMiniSat, CaDiCaL, Z3 and the pinned Z3 5.1.0 wheel
binary were obtained for items 3 and 5.

## Item 1 — XOR-learning CDCL (Xorricane / CDXCL) on the frozen `S′₄` corpus

- **Needs.** The Xorricane source (Danner–Kreuzer, JAIR 2026) and the upstream
  benchmark generator `mtrimoska/EC-Index-Calculus-Benchmarks`, or the
  repository's own regeneration of the 60-input corpus from
  `semaev_corpus.rs` (parameters are vendored; the ANF writer in
  `semaev_sat.rs` emits the system).  An ANF → XNF converter (Rust) is the
  one piece to write: each parity row with monomials becomes an XOR of
  AND-literals, which is 2-XNF directly.
- **Frozen inputs.** The 60 inputs of
  `research/index_calculus_baseline_20260914/regression/` at
  `(n, ℓ) ∈ {(15, 5), (17, 6), (19, 6)}`, plus uniform targets at
  `ℓ = 8, 9, 10` on `y² + xy = x³ + x² + 1` with `n = 3ℓ − 1`, 16 per rung,
  seed 20261011; `m = 3`.
- **Unit.** Conflicts per uniform target against `2^{3ℓ}/3!` (WDSat's
  attained bound, thesis Table 7.6) and against `T_enum = C(2^ℓ, 2)`.
- **Falsifier of the survey's verdict.** `< 0.5 · T_enum` at two consecutive
  `ℓ ≥ 8`.  **Abandon** if `≥ 0.5 · 2^{3ℓ}/3!` at `ℓ = 9` (no pruning above
  the last level).  Expected: abandon; class engineering at best.

## Item 2 — WDSat `-x` completeness audit, then SATIC

- **Needs.** WDSat at commit `61c6ff3` (the frozen regression's binary) and
  `55d55b2` (Stage 193), built by
  `research/index_calculus_baseline_20260914/pilot/build_pilot.py`; that
  runner is Python and must be replaced by a native driver before new
  measurements (`AGENTS.md`).  SATIC (Blomme–Cherif–Ionica–Dequen, CoDIT
  2025) if its code is published; the paper's full text was unreachable.
- **Frozen inputs.** The `solver_16` cell that returned UNSAT on a
  decomposable target at `n = 131` (`research/nagao_relations/solver_16/`),
  reduced to the smallest `d` that reproduces it under the 23 M-entry
  static cap; then the 60-input regression for SATIC.
- **Unit.** Verdict agreement with the exhaustive oracle (audit); summed
  median conflicts per configuration (SATIC vs WDSat symmetry-on, the
  regression's own target).
- **Falsifier.** Audit: the wrong UNSAT reproduces under `-x` and not in
  plain mode on the same instance → report upstream with the instance.
  SATIC: `< 0.5×` WDSat's summed median conflicts with no input `> 1.1×`
  → engineering; otherwise closed.

## Item 4 — Bosphorus ANF preprocessing before CryptoMiniSat and Kissat

- **Needs.** Bosphorus (Choo–Soos–Chai–Meel, DATE 2019) source; Kissat is
  not in the Ubuntu 24.04 archive (CaDiCaL 1.7.3 is).
- **Frozen inputs.** The same 60-input corpus in ANF; CryptoMiniSat 5.11.15
  and CaDiCaL 1.7.3 as the back ends, with and without `bosphorus --anfread
  … --cnfwrite …` (XL, ElimLin and learnt-fact simplification on).
- **Unit.** Conflicts on the 29 certified-UNSAT inputs, before/after.
- **Falsifier.** `< 0.5×` median conflicts after preprocessing.  Expected:
  no gain beyond the trace fall and the degree-2 rows the native path
  already injects; class accounting.

## Item 7 — GPU pairs-and-solve leaf

- **Needs.** A CUDA device for `gpu/fes/` (kernels compiled and
  cross-checked, never measured) or a Metal device.
- **Reference arm, measured here.**  The CPU pairs-and-solve oracle
  (`examples/semaev_decomp_bench`, one-word field, `gcd(q, L_V mod q)` per
  pair) on refuting targets to `ℓ = 11`:
  [`pairs_and_solve_cpu_l11.txt`](pairs_and_solve_cpu_l11.txt).  This is the
  arm a GPU kernel must beat.
- **Falsifier.** GPU candidates per second per device `< 50×` one CPU core
  at `ℓ = 12…16`.  Expected: engineering; it moves the toy-ladder reach,
  never the ratio (relation collection is `Θ(2^n)` with any enumerator).
