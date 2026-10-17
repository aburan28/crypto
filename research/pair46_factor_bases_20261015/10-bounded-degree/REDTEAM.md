# 10-bounded-degree — red-team review (2026-10-15)

Verdict: **not refuted.** The headline `exact_obstruction` is supported within its stated scope.

## Checks done
- **Preregistration.** PROTOCOL.md (mtime 16:27:13) is older than Cargo.toml, src, raw.json (16:30:57) and RESULTS.md. Its sha256 `e2b18e84...` matches PROTOCOL.sha256, raw.json and results.json. run.log shows a single run (87.6 s). The deviations are disclosed (RF run per prime, EXT checked by Hasse only, P46 seed mixing) and none of them matter.
- **Reproduction.** I rebuilt the instrument in a scratch copy and ran `taskset -c 3 ./target/release/bounded_degree 20000` (83 s). raw.json matched the committed file field for field. Only `wall_s` differed, and `protocol_sha256` was empty because the scratch copy had no PROTOCOL.sha256.
- **Independent recount (Python).** At p=65539 I got #E = 65736. For PW t^6 I got image 10924 and |F| = 10967, exponent 0.8388. All three match the instrument.
- **Numbers.** The results.json headline (min rho 0.167335, min exponent 0.838802, EXT max exponent 0.504977, gap 0.3338), the per-d c_d table, the claim that 13/14 Delta CIs contain 0 (RR d=7 [0.0002, 0.0023]), the P46 d=8 value 0.6411 [0.6282, 0.6540] and membership cost ≈157-160 d^2 all match RESULTS.md and the claim. The 168 rows are RP 63 + RR 63 + PW 21 + DK 21. The 63 ZD rows were checked separately (max |F| = 4 <= 2d).
- **Controls.** The identity map gives #E-1 exactly. The random function gives 0.6327 against 1-1/e. The x in F_q over F_{q^2} control gives |F| ≈ q, so the detector does see a sqrt-size factor base. PW and DK image counts are exact. The same Legendre-sum instrument is used on every arm.
- **Math.** The genus formula `br/2-1` reproduces floor((3d-1)/2) for deg h = 3d and 2d-1 for deg h = 4d-1. In X2, the map (t,Y) -> (P/Q, Y/Q^2) is onto F_f, with fibres <= d and the stated point losses, so LB_f is correct. In X3, (p-4d sqrt p)/d > sqrt p iff d < sqrt p/5, which gives 1.54e7 and 8.2e8 bits for pair 46. In X4, the pole order is max(2d, 2d-1) = 2d.

## Minor caveats (do not change the verdict)
1. **"Exponent" at toy sizes.** "exponent >= 0.839" and "gap >= 0.33" are log|F|/log p values for |F| ≈ p/d, that is, a constant fraction. They go to 1 as p grows, so at pair 46 the exponent would be ≈ 0.98. The quantity that carries over to larger p is rho (a constant ≈ c_d or 1/d). The exponent figures are a finite-p presentation of that constant, not a measured exponent. They are conservative, not inflated.
2. **X3's "d >~ p^0.1/5" is looser than needed.** The same bound gives d >~ p^0.1 (since p^0.9 >> 4 sqrt p). The weaker statement is still true.
3. **Seeds and families.** Pair 46 is one polynomial per d, sampled with n = 20000. The toy agreement with c_d is empirical (n = 9 per cell), and the 14 tests were not corrected for multiple comparisons. Both points are disclosed in RESULTS.md.
4. **Scope.** The claim covers only x-images of degree <= 8 maps and the A + yB conditions. This is stated correctly in RESULTS.md and must not be generalized.
