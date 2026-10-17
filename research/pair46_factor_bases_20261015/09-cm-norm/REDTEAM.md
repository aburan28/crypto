# 09-cm-norm: red-team review (2026-10-15)

**Outcome: the claim stands (not refuted).** The issues below are minor wording or scope points.

## Checks done
- **Preregistration.** PROTOCOL.md (mtime 16:17:56) is older than raw.json and run.log (16:19:47). Its sha256 is a1e2813a..., which matches PROTOCOL.sha256 and results.json. The arms, sizes, seeds, curves and decision rule were followed as written. Protocol C1 predicted v*lambda = 2 mod N. The run found this false, and RESULTS reports it openly as a correction.
- **Reproduce.** I rebuilt from a copy of src/ and ran with `taskset -c 3`: 15.9 s. The new raw.json is **byte-identical** to the committed one. run.log differs only by the `time` lines.
- **Numbers.** The ratios and CIs in results.json match RESULTS.md and the headline (m2 B256: 0.0159 [0.0156, 0.0161], AP 0.0082; m2 B1024: 0.0079 [0.0064, 0.0089], AP 0.0040; m3 B128: 0.0020, AP 0.0007).
- **Controls.** The null is uniform random symmetric sets of equal size on the same curve; this is the right null. The positive control (AP) passes its rule (≤0.1 at B≥256). The coverage function is the same on every arm.
- **Pair-46 math, checked independently in Python.** N = p+1-t, p = 5t^2+19t+19, 4p = t^2+19v^2, and p, N, v are prime. lambda^2-lambda+5 ≡ 0. v*lambda ≡ t and v*(1-lambda) ≡ 2. pi = (t-v)/2 + v*w = -1 + v*w. So pi acting as 1 forces the eigenvalue of w to be 1-lambda, as RESULTS says. Also 19v^2/4 = 5650313161883992 ≈ 0.95p, and h(Z[pi]) = v+1 (because (-19/v) = -1 and O_K* = ±1). The norm-form identity for a+b*pi is correct.
- **E1.** The sumset bound F_X+F_X ⊆ F_{4X} is correct. The measured |S+S| on curve 0 at B=256 is about 1991, close to 4|S| = 2048.

## Minor issues
1. what_was_cut says the toy lambda was "taken from the logs of the points". In the code, lambda is the two roots of x^2-x+5 mod N, computed from sqrt(-19). It is never tied to a real endomorphism. The ideal arm therefore depends only on N, not on the curve. This does no harm, since both roots run and give identical coverage.
2. The statement "ratio scales as ~4/B" holds only before the null saturates. At B=1024 on the 2^20 curves, null coverage_2 is about 0.86, and the measured 0.0079 is about 2x the value 4/B = 0.0039.
3. In E1, "decomposition costs Theta(sqrt N)" is a generic or heuristic statement, not an exact lower bound. The exact part is the sumset bound. A cleaner way to put the obstruction: the logs of the factor-base elements are known, so one relation solves the DLP. The chance that a trial succeeds is |F_{m^2 X}|/N, so this reduces to a BSGS-type time/memory tradeoff. The verdict still holds.
4. The bootstrap over curves mixes N ~ 2^20 and N ~ 2^22, so the B=1024 CI partly reflects that mix. It does not affect the decision.
5. The degree-5 endomorphism was not evaluated on pair-46 points (this is disclosed). The verdict does not depend on it.
