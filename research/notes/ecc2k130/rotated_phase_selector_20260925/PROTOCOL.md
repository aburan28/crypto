# Frozen n19 Frobenius-phase selector: exact support and oracle probes

Status: protocol and source freeze before any phase-closure or phase-query outcome. This is an exact archive analysis of the merged [n19 rotated PDP corpus](../rotated_pdp_corpus_20260925/RESULT.md), not an implicit S6 solver, relation-collection result, recovered logarithm, or speed comparison. The complete source histogram was independently replayed in #767; this experiment must independently replay its own phase transformation and selected row semantics. Do not use the archived target-by-target phase outcome to change either order.

## Distinct hypothesis

Use (E:y^2+xy=x^3+1) over (GF(2^{19})) with polynomial (x^{19}+x^5+x^2+x+1), prime subgroup order (q=130873), cofactor 4, (H=(385982,301867)), and (	au(H)=[41811]H). The factor product is the fixed beta-3 (m=6,d=2) construction (F_0\times\cdots\times F_5), one point from each distinct rotated slot. Global Frobenius sends this constrained product to a different six-slot window, so the older no-hit-rate-gain argument for a sigma-closed *admissible tuple family* does not apply. Its exact projected support (A) has 62,389 subgroup targets including O. Since (q-1=19\cdot6888), (lambda=41811\neq1) has order 19, and (62389\not\equiv1\pmod {19}), this (A) cannot be Frobenius-stable. That proves some added support, not a useful amount.

For target (Q=[j]H) and phase (s), query the **same** archived six-slot PDP at (	au^sQ=[\lambda^s j]H). The primary fixed order is (s_t=5t\bmod19), (0\le t<19); its spread is chosen structurally because neighboring slot windows share five of six spaces. The separately labelled control is contiguous (s_t=t). Each order has fixed caps (K=1,4,8,19); no order or cap may be selected after seeing target membership. A full 19-phase scan tests the same set in a different order, so both full-closure counts must agree.

For each order and every (j\in[0,q)), store the first successful **attempt position** as one byte (255 if none), its phase exponent from the frozen order, and all cap results. Define (C_K=\{Q:\tau^{s_t}Q\in A\text{ for some }t<K\}), with (C_0=\varnothing). Report (|C_K|), exact misses, first-hit-position counts, and the free perfect-oracle probe total
[
P_K=\sum_{t=0}^{K-1}(q-|C_t|)
   =\sum_{j=0}^{q-1}\min(K,\mathrm{firsthit}(j)+1),
]
where a full miss contributes (K). Count O and also report nonzero counts. These are complete finite-population counts, not a Bernoulli or solver-time model. The spread cap-4 comparator is #775's fixed four-beta union ((121851,249827)) for (support, four-probe total); the cap-19 context is its five-beta ((126265,258849)). Report the contiguous control and cap-8 results independently. Four/five beta arms have 12/15 potential signed columns; every phase arm retains the same three nonzero beta-3 columns.

The preregistered **support-and-probe priority screen** passes only if the primary spread order either (a) has cap-4 support at least 121,851 **and** at most 249,827 oracle probes, or (b) has cap-8 support at least 126,265 **and** at most ((9/4)q) oracle probes. The second bound is tested in exact integers as (4P_8\le9q). If both fail, retire this phase selector at n19. A passing screen only admits a separately frozen, fully charged complete-S6 solver comparison; it does not estimate ECC2K-130 attack speed. The contiguous order is a predeclared control, never a replacement for a failed primary order.

## Full-point and coefficient checks

Verify the group order, primality, all four rational 4-torsion points, (lambda^2+lambda+2=0\pmod q), exact order 19, and (	au(H)=[\lambda]H) with actual group arithmetic. Read the exact #767 beta-3 n19 projected histogram and six frozen factor lists. Check archive hashes, 117,649 total labelled tuples, 62,389 distinct projected points, and six seven-point factors satisfying (	au(F_i)=F_{i+1}). Build a complete (j\leftrightarrow[j]H) index; map every archived projected point (R=[4]Q) to exactly one (Q). The independent verifier uses separate bit-serial/Fermat arithmetic and rebuilds factor points directly from beta 3.

Use all 64 unconditional point-only holdouts frozen by #779. They were selected without consulting support and are not used to choose phases. For the first full-19 hit in each order, independently check the archived first labelled witness ((P_0,\ldots,P_5)): its full sum (S) equals (	au^sQ+T) for one exact (T\in E[4]), and ([4]S=[4]\tau^sQ). Undo slot (i)'s Frobenius on (P_i), project by 4, normalize a canonical signed column (C), and form a coefficient (epsilon_i\lambda^i). Check the aggregated row by group law,
[
\sum_i[\epsilon_i\lambda^i]C_i=[4]\tau^sQ,
\qquad
\sum_i[\epsilon_i\lambda^{i-s}]C_i=[4]Q.
]
Verify every observed (C_i) belongs to the same frozen three nonzero signed ([4]F_0) columns and that O contributes zero. Keep every selected hit/miss, phase, torsion index, witness index and row in the evidence. These checks do not solve the columns or recover a holdout log. The 64 holdouts are a fixed semantic sample, not an estimate of all-q support.

## Freeze, costs, and failure semantics

Before any outcome, commit this protocol, literal input, producer, independent verifier, runner, hash-only CI, and (FROZEN.json) pinning each source/input and the #767 archive, arithmetic references, and #779 point-only file. Run hash-only CI at the exact PR head and obtain peer review before one outcome run. The hash-only path may hash archives and compile source but must not extract projected histogram rows or calculate any phase membership. The producer has a 120-second wall and 512-MiB RSS cap; independent replay has 300 seconds and 512 MiB. Retain commands, UTC times, child wall/CPU/RSS, stdout/stderr, source/input/output SHA-256, binary first-hit arrays, all-q deterministic summary, fixed-case row certificates, and failed/censored partial evidence. CI after outcome must verify all hashes and independently recompute both complete first-hit arrays and sample rows. Do not repair an outcome by changing order, input, cap, or target list.

The possible charged successor must freeze a complete PDP exporter/solver with model lifting and externally checked negative proofs, then use the same 64 Qs and fixed phase schedule. Charge target Frobenius transforms, one-time base/log training, every phase export, SAT search, all completed UNSAT and unknown attempts, first verified witness, and independent point recovery. Compare same-Q four-beta and automorphism-aware batched rho at full end-to-end cost. An archive-oracle probe is never a measured solver query.
