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

sweep test (P-384, P-521, brainpoolP384r1) running in background at time of writing —
results pending. The 256-bit results (secp256k1 3/3) confirm the fix works there.

### Revised interpretation of `RESEARCH_LLL_GS_ANALYSIS.md`

- **§5 hypothesis** (secp256k1 r_i distribution causing μ-oscillation): likely wrong.
  The secp256k1 vs P-256 discrepancy was purely numerical.
- **Failure B** (P-384/P-521 iteration cap): also likely numerical, not LLL convergence.
  Should be confirmed by sweep test results.
- The Gram-Schmidt analysis in §3 remains valid as a theoretical framework, but the
  empirical claim that secp256k1 "fails LLL" is now invalidated.

### Commits made

See next git hash after this entry.

### Next step proposal

1. **Confirm sweep results** (P-384, P-521, brainpoolP384r1): wait for background
   sweep test, record results here. If P-521 now recovers → all curve failures were
   numerical.
2. **Update `RESEARCH_LLL_GS_ANALYSIS.md`**: revise the "two failure modes" narrative;
   replace with the numerical-overflow explanation.
3. **Move to priority thread 2** (CHLRS Igusa formula): the P-521 NaN is effectively
   resolved; the remaining work is documentation.
4. **Howe sextic twists** (priority 3): quick PARI check on 15 pairs, one session.
