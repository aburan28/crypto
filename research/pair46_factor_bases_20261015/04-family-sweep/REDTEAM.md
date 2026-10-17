# 04-family-sweep: red-team review

Reviewed 2026-10-15. The reviewer changed no other file in this directory.
Verdict: **the null verdict is supported, as scoped.** The issues below are caveats; none of them refutes a claim.

## Checks performed

- **Protocol freeze.** `sha256sum PROTOCOL.md src/main.rs` matches `PROTOCOL.sha256`
  (033bff61…, dae32149…). The freeze timestamp is 15:37:47Z, run.log starts at 15:37Z and its mtime is 15:50.
  main.rs (mtime 15:33) is the frozen version. The ordering rests on mtimes and the hash file, because there is
  no git commit. Arms, n, T, metrics, tests, α and the add-on all match the protocol.
- **Reproduction.** The full run exceeds the 10-min budget (749 s wall), so a partial rerun was done instead.
  The frozen main.rs was rebuilt in a scratch directory and run with `taskset -c 3 family_sweep out 120 20 26 99 50`,
  which took 57 s. Per-curve seeds depend only on the family index, so these rows can be compared directly:
  all 192 rows of `curves.csv` with bits ≤ 26 (F, N1, N2, PC) are byte-identical. The pair-46 add-on also
  reproduced exactly: crater 0.99881, ranks 73/200 and 82/200; floor 1.00250, rank 150/200. Startup checks
  passed: S3 self-test, pair-46 identity, j and order checks.
- **Numbers.** The arm means recomputed from `curves.csv` match the means in RESULTS.md and results.json
  (F 0.99275/0.62933/1.00076, N1 1.00873/…, N2 1.00779/…, PC 1.04854/…). M1 = decomps/(T·2|F_x|²/N)
  recomputed for all 600 rows matches the stored values. Every ratio, CI and p in the headline matches results.json.
- **Arms.** All 100 F rows satisfy p = 5t²+19t+19. None of the 200 N1 primes is in the family: we recomputed
  4p = a²+19y² for each and checked y ∉ {a+2, |a−2|}, and |t| = a. The N1 and N2 primes are all distinct.
  The same `measure()` function is used for every arm. Only PC passes dope = 1.03, which is the planted effect.
- **Math.** t² − 4(5t²+19t+19) = −19(t+2)² is correct. j(−19) = −96³ = −884736 is correct. The construction
  a = 3k, b = 2k with k = j/(1728−j) gives j (the code asserts it). Aut = {±1} holds because j ∉ {0, 1728}.
  Pair 46 has t = 34489701. The S3 polynomial in the code is the standard one, and the self-test exercises it.
  The normalizer 2|F_x|²/N is the right expectation for unordered x-pairs.

## Caveats (not refutations)

1. **Most of the M1/M2 positive control acts on the instrument, not the curve.** Widening the membership bound
   by 3% shows that the detector can see a ≈4% change in yield. It does not show that a curve-level structural
   effect would come through the instrument. The M3 control (best-of-2 model selection) is curve-level.
2. **The PC arm is not independent of N2.** Each PC curve is a re-model of the first N2 curve of the same
   member: same p and N, and isomorphic. The PC-vs-N2 permutation test treats the two arms as exchangeable,
   so the PC p-values are approximate. With p ≤ 0.0005, detection is not in doubt.
3. **Randomizing the model largely fixes the outcome in advance.** Each F curve gets a uniformly random model
   (a c⁴, b c⁶). The small-x factor base then depends on Legendre symbols of a random rescaling, so a near-null
   result is close to expected from character-sum heuristics. The run tests "random models of family curves".
   Only the M3 add-on tests the *specific* pair-46 models, and only on M3.
4. **Resolution.** "No family-level advantage" means no effect larger than about ±2% in M1, ±0.01 in M2 and
   ±0.0025 in M3 (the CI half-widths). RESULTS.md states this.
5. **Per-band results are underpowered.** The M1 positive control is *not* detected in the 25-27 band
   (ratio 1.0097, p = 0.64) or in the 20-24 band at α/6 (p = 0.010). A null per-band result therefore carries
   little information. RESULTS.md labels the per-band results descriptive, which is correct.
6. **The floor check is partial.** The run verifies #E = N and j ≠ −96³ for the given floor model. It does
   not verify End = Z[π]. That property comes from earlier items and is not established here.
7. **The disclosed deviation is handled correctly.** The expected PC M1 ratio is 1.061/1.012 ≈ 1.048,
   and the observed 1.040 [1.020, 1.061] agrees with it.

Nothing is extrapolated to pair-46 scale beyond the M3 add-on. No constant is presented as an exponent.
