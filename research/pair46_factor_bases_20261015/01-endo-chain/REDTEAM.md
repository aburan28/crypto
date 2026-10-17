# REDTEAM — 01-endo-chain (independent review, 2026-10-15)

Verdict on the claim: **not refuted.** The negative verdict holds. The decisive comparison
(endo-chain-known against bsgs-pm) uses two preregistered arms on the same instrument with
separated CIs, and exact identity E3 explains it.

## What was checked
- Protocol timing: PROTOCOL.md mtime 15:19:18, Cargo.toml 15:20, results.json 15:23. Its sha256
  `0ce98a09…dfcf2` equals `protocol_sha256` in results.json. All 4 sizes x 10 arms x 40 targets and
  phase 2 were run as specified.
- Re-run (phases 0-2, scratch copy of the binary, 21.8 s, 1 core): every count-based number reproduces.
  At 22 bits: bsgs-pm 0.955 [0.863,1.049]; nofold/k16 40.85; phase-2 rates identical to the 4th
  decimal. endo-chain-known came out 1.522 [1.320,1.759], against 1.47 reported. The only change
  is c_E (1.082 now), which is wall-clock timing, as the RESULTS file says it would be.
- E1 recomputed independently: lambda^2 - lambda + 5 = 0 mod N, and p = 5t^2+19t+19 with
  t = 34489701. t^2 - 4p = -19(t+2)^2. N-1 = 2*7*3037*139886589119, and ord(lambda) = (N-1)/7
  = 849671142308806.
- Instrument: every arm uses the same affine add, the same counters and the same verify_d (a
  charged scalar multiplication). The tables are built by the same chain_table, and the nofold
  positive control reuses the k16 table with only the decoder changed. This is the correct control.
- E3 cost model: the measured known/bsgs ratios (1.41-1.63) match sqrt(2 c_E) = 1.46-1.51.
- E4: j = -96^3 is not 0 or 1728, so Aut = {+-1}. Correct.
- E5: hits for u occur at a = u - d and a = -u - d. Correct, and Theta(N) for a sequential walk.

## Issues (none overturns the verdict)
1. **The rho-neg reference is defective, and RESULTS mis-describes why.** Every useless collision
   triggers a fresh start costing two scalar multiplications. With r = 20, a fruitless 2-cycle comes
   about every 2r steps, so this overhead is proportional to the walk length times log N; it is not a
   constant. At 22 bits it is about 45 restarts x about 67 ops, roughly 1.46 sqrt(N). Adding about
   0.65 sqrt(N) of set-up and 0.88 sqrt(N) of walk reproduces S = 2.93. The sentence "the set-up is a
   constant ..., so rho's S tends to its walk value" therefore holds for rho, but not for rho-neg.
2. **The primary effect size was swapped without being flagged.** The protocol's effect size is the
   ratio to rho-neg mean S, and its "advances" rule is stated against rho-neg S. Taken literally,
   endo-chain-known has its CI upper bound below rho-neg's CI lower bound at every size. "advances" is
   blocked only by the "ratio does not grow" clause: the ratio grows from 0.31 to 0.50, and that growth
   is an artefact of issue 1. RESULTS reports bsgs-pm and the walk-only rho-neg steps instead, and does
   not tabulate `ratio_to_rho_neg`. That substitution is sound, because BSGS is a preregistered and
   stronger reference on the same instrument. It should still be labelled a deviation, and the
   preregistered ratio should be shown.
3. **The protocol's E3 was silently corrected.** PROTOCOL.md says ordinary +-BSGS costs sqrt(2N).
   RESULTS says sqrt(N), which is the right average-case cost for m = sqrt(N)/2 and is what the
   measurement gives. Under the protocol's formula the predicted known/bsgs ratio would have been
   sqrt(c_E), about 1. The correction is right but is not noted.
4. **"About sqrt(k)" is loose, and "constant" holds only at fixed k.** In the post-hoc run the 22-bit
   ratios against the abscissa base are k4 0.42 (sqrt(k) predicts 0.50), k16 0.29 (0.25), k64 0.19
   (0.125) and kF 0.096 (0.022). At fixed k the ratio still drifts down with size: k16 goes 0.53,
   0.44, 0.34, 0.29. So the constant has not been reached by 22 bits, and only the heuristic makes the
   result a constant. If k grows with N (k = F), chains improve on the abscissa base by an exponent,
   from about N^{3/4} to N^{1/2}. They still do not beat BSGS (E3). The headline's "at fixed |FB|"
   wording is acceptable only when read as "fixed k" as well. These runs are post-hoc and are labelled
   as such.
5. **Phase 2 has limited power.** The arithmetic-progression positive control is an extreme
   structure (about 100x off the null). It shows the detector works, but not that it could see subtle
   psi-induced coverage loss. "No structure detected" covers m = 2, one curve and one |FB|, which is
   how RESULTS scopes it.
6. **There is no overstatement for pair 46.** No toy S value is extrapolated to pair 46. The pair-46
   claims are E1 and a measured c_E of 0.98 (not isolated, host shared, as disclosed). The c_E of
   about 2.2 under batched inversion is labelled as derived.
