# Full-rank column certificate: implementation gate, no new speed claim

The [n37 full-rank panel](../full_rank_compact_orbit_ic_20261003/RESULT.md)
reached rank 43/43 and verified every recovered target logarithm. Its matrix
therefore fixes all 42 factor-base unknowns, but the run did not separately
check their values against curve points. This PR adds the opt-in
`incremental-gauss-full-rank-checked` mode to close that correctness gap.
After the unchanged full-rank relation stream recovers the target, the new
mode checks one point with a nonzero coefficient from each column using
`[coefficient × column_log]G = [cofactor]P`. Both multiplications, misses,
and failures enter the verification ledger. A failed or missing column
prevents a verified result. The existing full-rank mode and its archived
session retain their original method identity and charge.

This is a certification feature, **not a new n37 measurement**. The focused
Rust tests compare the checked and historical modes on the same prime-curve
relation stream, reject a false column log, and reject a declared column
without a representative. The Linux `ecbench` CI gate replays the archived
n37 sessions. No IC/rho speed ratio or ECC2K-130 transfer is inferred from
these tests.

The evidence-ranked next step is a **shared-table, shared-rank, checked
batch** on new orbit-disjoint n37 public points against a same-point strong
automorphism-aware batched rho process. Its protocol must fix batch size,
rank and target seeds, native implementation, timeout and memory caps,
complete cold and per-target work, and an independent full-point replay
before either arm runs. Charge column certification once with the rank
database, every failed descent, and all target checks. The earlier
[1,024-target cold panel](../disjoint_cold_v2_outcome_20261001/RESULT.md)
already rules out repeating the old W64 policy; the newer folded table is
the specific unmeasured variable. Likewise, the
[equal-size n41 four-window screen](../base_window_screen_outcome_20261001/RESULT.md)
found no selector opportunity among those windows, so another same-window
screen is lower priority.

The [degree-263 ring certificate](../../../ecc2k130_endo_ring_263_20260925/README.md)
establishes a conductor-263 descending order with discriminant
`-7 × 263² = -484183`, but no cheap squaring endomorphism on one fixed leaf.
That supports a separately costed descendant-native high-arity PDP gate;
it does not predict yield or let the leaf inherit the source's free orbit
action. Changes to SAT, F4/F5, FES, crossbred or Gray-code kernels should be
ranked by complete PDP and recovery costs after this gate, rather than by
isolated solver timings. The n37 four-policy PDP sample already saturated
at three summands for every tested policy, leaving no yield selector signal
there; larger degrees and the degree-263 leaf remain open.
