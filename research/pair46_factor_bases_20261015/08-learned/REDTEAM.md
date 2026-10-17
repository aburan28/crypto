# 08-learned — red-team review (2026-10-15)

Verdict: **claim supported, with minor overstatements.** Not refuted.

## Checks performed
- Protocol timing: PROTOCOL.md mtime 16:11:18; run_learned.py 16:12:18; results.json 16:14:37. sha256 of PROTOCOL.md matches PROTOCOL.sha256 and results.json `protocol_sha256`. The design (bits, curves, seeds, mu, T=N/16, arms, features, L1 settings, bootstrap) matches the code. The SMOKE switch only restricts BITS.
- Reproduction: `SMOKE=1 OMP_NUM_THREADS=1 python3 run_learned.py <scratch>` (14-bit, 24 units, 4 s) gave output **byte-identical** to the 24 14-bit units in results.json.
- Instrument: every arm uses the same `counts_all` FFT count, the same training mask and the same held-out set (r=0 excluded). Bases are k in 1..(N-1)/2, which is equivalent to choosing x-values. Random arms are uniform over x.
- Numbers: summary.json reproduces the RESULTS.md table, the in-sample figures (greedy 0.92/0.82/0.76 vs random 0.64 at mu=1), the greedy held-out mean-count ratio 0.91-0.98, PC1 mean-count ~1.0, the AUC cell means (0.48-0.55 and 0.50-0.51), PC2=1.000, and the identity check on 504 arm evaluations.
- Math: sum_R c(R)=|S|^m is trivially exact for any S in any group, so E[c]=|S|^m/N for every base. |S+S| <= L(2L+1) for a symmetric S. Random success matches 1-e^{-mu} (0.221 and 0.632 observed). The ceiling mu/(1-e^{-mu}) gives 1.13 and 1.58.
- Controls: PC1 and PC2 pass. random2 excludes 0 at 16 bits at both mu. This is reported honestly as making the empirical part formally inconclusive under the frozen rule. Pooled over 72 units, random2 vs random is +0.0009 +- 0.0015, which is consistent with the miss being an interval-coverage issue.

## Issues (minor)
1. Headline "greedy is 0.4-4.6% below random on held-out" is false in one cell: at mu=1, 14 bits greedy is +0.0006 [-0.005,+0.006]. Likewise "learners beat greedy (+0.005 to +0.051)" fails in that cell: learned_greedy +0.000 [-0.007,+0.006], learned_success +0.004 [-0.002,+0.010].
2. "Held-out success is at most L(2L+1)/N" is slightly wrong. A base can place its sums in the held-out set, so the held-out bound is L(2L+1)/(N-1-T), about 16/15 larger. The ceiling line in fig_heldout uses /N.
3. The headline wording "caps any selector's success advantage at L(2L+1)/N" conflates the success cap with the ratio cap. It is the success that is capped at L(2L+1)/N, and the ratio vs random is about mu/(1-e^{-mu}).
4. Calling mu -> 0 "the index-calculus operating point" is debatable: m-term IC usually runs at a constant per-target success, about 1/m!. The conclusion "constant-only ceiling, no exponent change" still holds. For m >= 3 only the mean identity is exact; no success ceiling is derived for m >= 3.
5. Scope: the obstruction is about relation *counts and success probability* for uniform targets. It says nothing about bases chosen to make the decomposition *solve* cheaper (structured or subspace bases for summation polynomials). So "exact_obstruction" should not be read as "factor-base choice cannot help IC on pair 46".
6. The per-unit learner AUC test sets are tiny (20% of 2L samples, e.g. about 16 at 14 bits), and per-unit aucG ranges 0.32-0.94. Cell means are fine, but per-unit AUC is noisy.
7. learned_success labels use `part > median` (hi) against `<= median` (lo), while the protocol says "at least the median". This is immaterial.
