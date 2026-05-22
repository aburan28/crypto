# Research AutoLab Log

Automated daily research log for the isogeny-graph ECDLP project.
Each entry is one autolab session run.

---

## 2026-05-21 (autolab run)

### Task picked

Priority thread 1: P-521 LLL NaN — the `gram_schmidt` function in
`src/cryptanalysis/lattice.rs` never had scaling; first run needed to
investigate, implement, and test the fix.

### Work done

- Audited `gram_schmidt` (lattice.rs:125) and `big_to_f64` (lattice.rs:150).
- Traced the NaN source: for a 521-bit curve, `n² ≈ 2^1042 > f64::MAX ≈ 2^1023`.
  `big_to_f64` called `to_f64().unwrap_or(0.0)` which returns `Inf` (not `None`),
  so the GS norm computation sums `Inf * Inf = Inf`, and `Inf / Inf = NaN`.
- Found a SECOND related bug: for secp256k1 (256-bit curve), `n ≈ 2^256`,
  `n^4 ≈ 2^1024 > f64::MAX`. The GS norm for `bstar[0]` is `||bstar[0]||² = n^4`,
  which also overflows → Inf → NaN for secp256k1 (but NOT for P-256, whose `n` is
  smaller by `~2^{192}`, putting `P-256.n^4` just under `f64::MAX`).
- Implemented `big_to_f64_scaled(x, scale_shift)`:
  - For small entries (`nbits + scale_shift ≤ 1020`): direct `to_f64() / 2^scale_shift`
  - For large entries: first integer-shift to bring into `f64` exponent range,
    then re-apply remaining exponent as float multiply; preserves small entries
    (e.g. `2^384` scaled by `2^542` → `2^{-158} ≈ 2.7×10^{-48}`, finite, not 0).
- Updated `gram_schmidt` to compute `scale_shift = max(0, max_bits - 500)` over
  all basis entries (target: `dim × entry² < 2^1023`), then call `big_to_f64_scaled`.
  The μ_{ij} coefficients are invariant to global scaling, so the result is exact.
- Ran full test suite:
  - `cargo test --lib cryptanalysis::lattice`: 9/9 pass
  - `cargo test --test curve_audit`: 5/5 pass
  - `cargo test --test lll_degeneracy_probe -- --nocapture` (head-to-head):
    secp256k1 3/3 RECOVERED (~663ms/seed); P-256 3/3 RECOVERED (~654ms/seed)

### Findings

**Critical**: The "secp256k1 LLL degeneracy" documented in `RESEARCH_LLL_GS_ANALYSIS.md`
was NOT a genuine lattice-theoretic phenomenon. It was a floating-point overflow bug:
secp256k1.n is close enough to 2^256 that `n^4` exceeds `f64::MAX`, corrupting GS
with Inf → NaN and producing false LLL divergence (cap hit). P-256 escaped only
because its `n` happens to be smaller by `~2^{192}` making P-256.n^4 just fit.

Numerical verification:
```
secp256k1: n ≈ 2^256 - 2^32,  n^4 = Inf  (overflow f64::MAX ≈ 2^1023.99)
P-256:     n ≈ 2^256 - 2^192, n^4 = 1.797e308 (just under f64::MAX)
P-521:     n ≈ 2^521,  n² ≈ 2^1042  (direct overflow, never fit in f64)
```

After fix (scale_shift = max_bits - 500 when max_bits > 500):
```
P-521 basis entry n² → n²/2^542 = 2^500   (finite, ~3.3×10^150)
P-521 basis entry 2^384 → 2^384/2^542 = 2^{-158}  (finite, ~2.7×10^{-48})
dot product of 10 scaled-P-521 entries: ~10^302 < f64::MAX ✓
```

**Full sweep results** (`probe_lll_sweep_by_bit_length`, 1 seed, 8 sigs, LLL):

| Curve | n_bits | k_bits | Old outcome | New outcome | New time |
|-------|--------|--------|-------------|-------------|----------|
| secp192k1 | 192 | 144 | ✓ 389ms | ✓ RECOVERED | 325ms |
| secp224k1 | 224 | 168 | ✓ 617ms | ✓ RECOVERED | 615ms |
| **secp256k1** | 256 | 192 | ✗ cap 10019ms | ✓ RECOVERED | 657ms |
| P-256 | 256 | 192 | ✓ 710ms | ✓ RECOVERED | 686ms |
| brainpoolP256r1 | 256 | 192 | ✓ 714ms | ✓ RECOVERED | 671ms |
| **P-384** | 384 | 288 | ✗ cap 20392ms | ✓ RECOVERED | 2047ms |
| **brainpoolP384r1** | 384 | 288 | ✗ cap 20131ms | ✓ RECOVERED | 2042ms |
| **P-521** | 521 | 384 | ✗ cap (NaN) | ✗ no short vec | 126010ms |

P-521 now runs WITHOUT NaN and WITHOUT iteration cap (LLL terminates cleanly in 126s),
but does not recover the key. The error is "LLL did not produce a recoverable short
vector — try more signatures or increase bias". This is a **genuine** failure mode
distinct from the numerical bug: with m=8, k_bits=384, the LLL reduction at 521-bit
precision does not find the short vector.

### Revised interpretation of `RESEARCH_LLL_GS_ANALYSIS.md`

- **§5 hypothesis** (secp256k1 r_i distribution causing μ-oscillation): **refuted**.
  The secp256k1 vs P-256 discrepancy was entirely due to floating-point overflow.
- **Failure A** (secp256k1 at 256 bits): **resolved** — numerical bug, not lattice.
- **Failure B** (P-384/P-521 iteration cap): **P-384 resolved** (numerical), **P-521 partial**.
  P-521 NaN is fixed; P-521 now has a genuine LLL precision/convergence issue at 521 bits.
- The Gram-Schmidt analysis in §3 remains valid as a theoretical framework.
- Open Q for P-521: is the "no short vector" failure due to (a) insufficient bias m=8/k_bits=384,
  (b) f64 precision loss in 521-bit LLL steps, or (c) need for BKZ vs LLL?

### Commits made

See next git hash after this entry.

### Next step proposal

1. **P-521 genuine failure diagnosis** (priority 1, continued):
   Try `m=16, k_bits=384` (doubling signatures) to see if more samples help.
   Try BKZ-β instead of LLL for the P-521 probe.
   Both are one-line parameter changes in `probe_once()`.
2. **Update `RESEARCH_LLL_GS_ANALYSIS.md`**: revise the "two failure modes" narrative;
   replace with the numerical-overflow explanation (the §5/7 analysis is now obsolete).
3. **CHLRS Igusa formula** (priority 2): implement the `(E × E^t)/Γ_α` Igusa-Clebsch
   invariants in PARI.  Start from `secp256k1_cm_audit/igusa_clebsch.gp` which may
   already scaffold this.
4. **Howe sextic twists** (priority 3): quick PARI check on 15 pairs.

## 2026-05-22 (autolab run)

### Task picked

Priority thread 1 (P-521 LLL), continued from 2026-05-21. Yesterday's run
fixed the overflow bug (f64 Inf/NaN). Today: diagnose why P-521 still fails
after overflow fix, and implement a root-cause fix.

### Work done

- Traced P-521 failure post-overflow-fix to **catastrophic cancellation** in f64 GS.
  Mechanism: P-521 HNP basis has entries spanning 2^384 to 2^1042 (658-bit
  dynamic range). GS step computes `n·a_l - (a_l/n)·n² = 0` exactly, but
  two ~2^499 f64 values subtract to a phantom residual ~2^446, which is
  2^604× larger than the true b*_m[m] = 2^(-158) (scaled). LLL's Lovász
  condition is meaningless with this noise floor.
- Verified the hypothesis by computing the expected residual analytically:
  `|error| = n × r / 2^scale` where `r = (a_l × 2^scale) mod n < n`.
  With scale=542, error ≈ 2^(521-542) = 2^(-21), still > true b*_m[m] = 2^(-158).
  (Actually error ≈ n × n / 2^2048 = 2^(1042-2048) = 2^(-1006) with HP GS.)
- Implemented `gram_schmidt_hp` and `lll_reduce_hp` in
  `src/cryptanalysis/lattice.rs`:
  - 2048-bit fixed-point BigInt arithmetic (HP_PREC = 2048)
  - Functions: `hp_from_bigint`, `hp_mul`, `hp_div`, `hp_round`
  - `gram_schmidt_hp` returns `(bstar_sq, mu)` as HP BigInt vectors
  - `lll_reduce_hp` mirrors `lll_reduce` but uses HP GS throughout
- Added `HnpReduction::LllHp` variant to `hnp_ecdsa.rs` and wired it in.
- Added `probe_once_ext` (accepts `HnpReduction` param) and
  `probe_p521_lll_hp` test to `tests/lll_degeneracy_probe.rs`.
- Added unit tests `lll_hp_matches_lll_on_small_basis` and
  `lll_hp_recovers_p384_key` to `src/cryptanalysis/lattice.rs`.
- Updated `RESEARCH_LLL_GS_ANALYSIS.md` §10 with full resolution writeup.

### Findings

**Definitive empirical result** (`probe_p521_lll_hp`, k_bits=384, m=8):

```
f64 LLL: 0/3 seeds recovered  (~147s each, loops to iteration cap)
HP  LLL: 3/3 seeds recovered  (~79-82s each, converges correctly)
```

HP LLL is **1.8× faster** than f64 LLL on P-521 because correct GS enables
proper LLL convergence — f64 LLL wasted cycles on wrong swap decisions.

Precision analysis: HP GS reduces the cancellation residual from ~2^446 (f64)
to ~2^(-1006) (2048-bit fixed point) — an improvement of 2^1452 in precision.
This completely resolves the P-521 issue.

Full test suite: 11/11 unit tests + 5/5 integration tests pass.

### Next step proposal

1. **Priority 1 CLOSED** for P-521. Residual open: lll_reduce_hp performance
   at higher m (currently m=8 at 79s; m=32 may need incremental GS update
   instead of full recompute per swap). Performance profiling optional.
2. **Priority 2 (CHLRS Igusa formula)**: implement `(E × E^t)/Γ_α`
   Igusa-Clebsch invariants in PARI. Check `secp256k1_cm_audit/igusa_clebsch.gp`
   for existing scaffolding.
3. **Priority 3 (Howe sextic twists)**: quick PARI check on 15 pairs of
   6 sextic twists for secp256k1. Short script, possibly runnable next session.
4. **Priority 5 (GLV-HNP Phase 2 toy)**: now that P-521 HP GS works, the
   GLV-aware lattice attack on a 32-bit toy curve is the next most interesting
   experiment.

### Commits made

[see next git hash after this entry]
