# l=8 subspace S₄ FFD measurement — binary decomposition frontier push #2

**Date:** 2026-09-20 · **Host:** Adams-MacBook-Pro.local (Apple M4 Pro, 48 GB, macOS 26.6, Darwin 25.6.0 arm64)
**Ledger gate served:** `docs/ic/boundary_targets.json` → binary decomposition `next_target`:
sub-2^(2ℓ) oracle at ℓ=8 — "ffd min/max/mean over >=16 draws present (or explicit inapplicable+reason)".
Baseline arm reference: `runs_manual/binary_decomp_l8_20260911/` (pairs-and-solve, 68.0 ms/target ledger-recorded; "FFD unmeasured" for the subspace system).

## What was built

- `src/cryptanalysis/ffd_harness.rs`: degree-general Macaulay measurement over monomial-set
  equations (`build_macaulay_rows_monomial`, `measure_monomial_system`,
  `macaulay_memory_estimate_monomial`), same conventions as the quadratic path
  (multilinear basis, truncation at D, multipliers ≤ D−δ per equation, operational
  fall = smallest D with rank < rows AND rank < cols).
- `examples/ffd_s4_subspace.rs`: builds the Weil-descended symmetrised S₄ system
  (`binary_semaev_s4::weil_descend_s4`, Koblitz b=1, n=24, ℓ=8, same irreducible
  `z²⁴+z⁴+z³+z+1` as the baseline arm), **exactly eliminates the 45 e-variables**
  (correspondence ANFs are target-independent — asserted at runtime), giving 24
  equations of degree ≤6 in the 24 factor-base bits, and measures two variants:
  - **full**: 24 eqs / 24 x-vars, D ∈ [4,8] (D=9 dense estimate ≈ 18 GB, skipped under 2 GB cap);
  - **x1-fixed**: X₁ frozen to a random V-element per draw → 24 eqs / 16 vars, degree ≤4,
    D ∈ [4,7] (D ≥ 8 is Koszul-ambiguous for δ_min=4: f_i·f_j − f_j·f_i enters at degree 8).

## Correctness gates (all passed, fail-closed)

1. **Sanity**: degree-general builder reproduces the quadratic harness rank-for-rank on
   descended S₃ at n∈{5,6,7} (seed 0xFFDDEAD), including documented FFD=3.
2. **Elimination exactness**: on 64 random V³ points per draw (1024 total), the eliminated
   system is satisfied **iff** `symmetrised_s4_eval(σ(X), x_R)=0` — zero mismatches.
3. **Witness vanishing**: both decomposable draws' witnesses (from `semaev_decomp::decompose`,
   independently attested) satisfy the eliminated system and its X₁-fixed fold at the witness X₁.

## Results (16 draws, seed 0x54FFD518, 2748.8 s total)

| variant | eqs/vars | δ | degrees measured | fall | notes |
|---|---|---|---|---|---|
| x1-fixed | 24/16 | 4 | 4–7 | **D=7 on 16/16 draws** (min=max=mean=7) | rank 15940 < rows 16728 at D=7; **deficit 788**; 7 < 2δ_min=8 → Koszul-clean ("structural") |
| full | 24/24 | 6 | 4–8 | none (rank=rows at D=6,7,8: 24/600/7224) | censored at D=9 (≈18 GB dense) |

- **Deficit invariance:** 788 on every draw — identical on the 2 decomposable and the 14
  non-decomposable targets ⇒ the D=7 syzygy space is **target-independent geometry**
  (expected origin: relations induced by the 16-dim σ-image of (X₂,X₃) inside the 37-dim
  (e₂,e₃) space), not a decomposability signal. It cannot decide the decomposition
  question directly.
- **Pairs recheck on this host:** median 16.5 ms/target (warm, this binary). The
  ledger-recorded baseline is 68.0 ms from the 20260911 arm; runs are reported side by
  side and never mixed.

## Obstruction (quantified, scope: ℓ=8, dense F₂ Macaulay LA, this host)

- Per-X₁ Gröbner oracle shape: Macaulay at the fall degree D=7 costs ≈30 s per X₁
  (rows 16728 × cols 26333, build-dominated) ⇒ ≈2.1 h per target over the 2⁸ sweep.
- Full-system shape: D=8 LA costs ≈135 s per target (7224×1271626 bits ≈ 1.15 GB) and
  shows **no** fall; column saturation (refutation certificates) needs rows ≥ cols,
  unreachable below ≈18 GB (D=9) and ≈760 GB (D=10).
- Baseline pairs-and-solve: 68 ms (ledger) ⇒ the direct Macaulay/XL route is
  **10³–10⁵× too slow** at this rung. The sub-2^(2ℓ) oracle gate is **unmet**; this run
  closes only the FFD-logging sub-gate.
- **Resource check (obstruction as asset):** the 788-dim, target-independent syzygy
  module at D=7 is a compressed description of the image geometry. A mechanism that
  quotients by it symbolically (F₅-style preprocessing, or precomputed syzygy reuse
  across X₁ — the module is target-independent, so it is computed once per curve/ℓ,
  not per target) could cut the per-X₁ row count materially; whether the residual
  system then decides below the pairs cost is open and testable at this rung.

## Forward guidance

1. Extract the D=7 syzygy module (kernel of the x1-fixed D=7 Macaulay transpose) and
   identify its generators' origin; verify target-independence directly (expected).
2. Test syzygy-quotient Macaulay (precompute once, reuse across X₁ and targets) — the
   only Macaulay-shaped route whose per-target cost could plausibly approach the gate.
3. Otherwise the beat's sub-2^(2ℓ) oracle stays open to non-Macaulay mechanisms
   (batched resultant/gcd structure over X₂, Frobenius-invariant factor bases —
   constant-factor, changes the comparison object — or higher summation shapes).

## Files

- `ffd_s4_subspace_summary.json` — full per-draw measurements (fixture_hash in claim draft pins it)
- `ffd_s4_subspace.stdout.txt` / `.stderr.txt` (empty) — run log
- `claim_draft_s4_subspace_ffd.json` — schema-v2 claim draft (measurement-only non-claims listed)
- `claim_check_s4_subspace_ffd.json` — autolab claim-check **PASS** (fail-closed)
- `make_claim_draft.py` — regenerates the draft from the artifacts
- `source_and_binary_hashes.txt`, `host.txt`, `launch_started_at.txt` — provenance

**Replay:** `cargo build --release --example ffd_s4_subspace && ./target/release/examples/ffd_s4_subspace --draws 16 --seed 0x54FFD518 --dmax-full 8 --dmax-fix 7 --mem-cap-gb 2.0 --out <fresh>.json` (~46 min). Deterministic given the seed; sanity gate aborts on any convention drift.

**Status:** PENDING_INDEPENDENT_VALIDATION. No ledger promotion is proposed: the oracle
gate is unmet and the measurement has not been independently replayed. No key-recovery,
sub-2^(2ℓ), or asymptotic claim is made; scope is ℓ=8 (n=24) only.
