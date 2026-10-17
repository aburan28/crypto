# Protocol (frozen before data): Joux–Vitse `(n − 1)`-decomposition at `n = 5`

Registered 2026-10-15, before any claim run.  The source state it freezes is
the commit that adds this file (`src/cryptanalysis/joux_vitse.rs`,
`examples/joux_vitse_ladder.rs`).  Runs before this commit were exploratory
probes (p = 11, 31, 241; a 16-trial DLP smoke test at p = 11) used only to
size the budget; none of their numbers is reported as a result.

## Question

Does Joux–Vitse's `(n − 1)`-decomposition (decompose `R ∈ E(F_{q⁵})` into
four factor-base points, `F = {P : x(P) ∈ F_q}`, through the overdetermined
descended `S₅` system) give a route where measured index-calculus cost
approaches or beats Pollard rho at sizes this repository can run?  If not,
how far off, and why?

## Boundaries (stated before measuring)

- **Reference: Pollard rho on the same group, same unit.**  The unit is `F_p`
  multiplications (Fp5's counter for curve arithmetic, f4_fp's counter for row
  reductions).  For rho, `fp_muls_per_op` is measured on the same `Curve5`
  arithmetic; the expected rho cost on a group of prime order `ℓ` is
  `√(πℓ/2)` group operations (r-adding walk, no negation map) or
  `√(πℓ/4)` with the negation map.  In the DLP arms rho is run on the
  same target and its count is measured.
- **Floor (counting): trials per relation.**  With `|F|` base points modulo
  `±`, at most `16 · C(|F| + 3, 4)` points of `E` are a signed sum of four
  base points (multisets, all sign patterns), so the decomposition rate is
  at most `ρ_max = 16 · C(|F| + 3, 4) / #E ≈ 2/(3 q)` for `|F| ≈ q/2`.
  Expected trials per relation are at least `1/ρ_max`.  The thesis heuristic
  is `1/ρ = (n − 1)! q = 24 q` (Vitse 2011, §7.3.2, Hypothesis 7.3.2).
  The measured rate is reported as a ratio to `ρ_max` and to `1/(24q)`.

## Fixed method (not tuned during the study)

- Curve: `y² = x³ + ax + b` over `F_{p⁵} = F_p[t]/(t⁵ − c)`, `p ≡ 1 (mod 5)`,
  `a, b` uniform with `j ∉ F_p`, drawn from `StdRng::seed_from_u64(seed ^
  (p << 20))` (ladder) or `seed ^ (p << 20) ^ 0xd1` (DLP).
- System: `H(e₁, …, e₄; x_R)` by interpolation (495 monomials × 9 powers),
  checked on 8 fresh points; five `F_p`-components.
- Trial: `f4_fp` F4, grevlex, degree bound 36 (Macaulay bound, Vitse's
  Hypothesis 7.3.4); trace recorded on the first inconsistent system of the
  curve and replayed (F4Remake) with fallback to full F4 on mismatch (both
  charged).  A consistent system is solved with `f4_fp::solve` (same options),
  split, lifted, and every sign pattern tried; a relation is kept only if it
  re-adds to `R` in the group.
- Oracle: meet in the middle over pairs (`mitm_decompositions`), no
  polynomial system.

## Runs

1. **Ladder** (`trials` mode), `p ∈ {11, 31, 41, 61, 71, 101, 131, 151, 181,
   191, 211, 241}`, seeds `{1, 2}`:
   - `--f4-trials 32`: sequential F4 trials, each cross-checked against MITM;
   - `--mitm-trials 500·p`: MITM-only trials for the decomposition rate.
   - Output `data/ladder/p<p>_s<seed>.jsonl`.
2. **Complete DLP** (`dlp` mode), paired with rho on the same `(P, Q)`:
   `p = 11`, seeds `{1, 2, 3}`, `--min-ell-bits 10`; `p = 31`, seed `1`,
   `--min-ell-bits 16`; `--max-trials 200000`.  Output
   `data/dlp/p<p>_s<seed>.jsonl`.  If the `p = 31` arm does not finish within
   10 hours of wall time it is stopped and reported as incomplete with its
   partial record.

Host: the session container (4 vCPU Intel Xeon @ 2.10 GHz, 15 GB),
isolation none; wall times are practicality notes only.  Counts are the
metric.

## Analysis plan

- `C_trial(p)`: mean over the 64 F4 trials of `muls_system + muls_f4 +
  muls_solve`, with min/max and sd; separately the replay-only and full-F4
  counts.  Fit `log C_trial` against `log p` (prediction: slope ≈ 0).
- `ρ̂(p)`: (MITM hit trials + F4 trials with a relation) / all trials,
  Wilson 95 % interval; MITM counts decompositions into four distinct base
  points only, so the F4 trials also report repeated-point relations.
- Cost per relation (derived, not measured as one quantity):
  `C_rel(p) = C_trial(p) / ρ̂(p)`.  Fit against `p` over the rungs with at
  least 5 hits (prediction: slope ≈ 1).
- Model (extrapolation, labelled as such): total IC cost
  `≈ (|F| + 1) · C_rel + LA`, with `LA ≈ 4 |F|²` multiplications modulo `ℓ`
  (sparse, 4 entries per row), against rho `√(πℓ/4) · c_add` with
  `ℓ ≈ p⁵`; solve for the crossover `q*` and compare with Vitse's stated
  `2^{60}` (§7.3.3).
- Comparison with the `n = 4` D₄ oracle per relation (PR #1863: `C₄ =
  5.5–8.5·10⁷`, about `4!` trials per relation), labelled as a comparison
  across different groups (`p⁴` against `p⁵`).

## Success and abandonment criteria (declared in advance)

- **Route found** if, in any DLP arm, `IC online F_p mults / rho F_p mults ≤
  10`, or if the measured `C_rel` and rate put the model crossover below
  `2^{80}` for `#E`.  Expected (from the thesis figures): neither holds;
  `IC/rho ≈ 10⁵–10⁷` at runnable `p`, crossover near `q ≈ 2^{60}`.
- **Instrument fails** if any F4 trial's distinct-decomposition set differs
  from MITM's, or any relation fails to re-add, or a DLP answer fails `kP = Q`.
  Such a run is reported, not dropped.
- Inadmissible: changing the degree bound, the trace policy, the factor base
  or the seeds after data; dropping setup or failed trials from the
  accounting.
