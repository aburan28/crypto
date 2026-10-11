# n=83 F6 wide-pair residual batch-size experiment

Registered before measurement on 2026-10-04. The hypothesis is that the
256-element residual batch in `F6WidePairIndex::solve4` pays for too many
field inversions and temporary vectors. Change only that batch size to 2,048;
keep pair construction, lookup, verification, and the public target identical.

Use the registered K0 curve `icv1-f2m83-tm6151469093347-debefd74`, its
standard polynomial-subspace cofactor-projected factor bases at dimensions
8, 10, and 12 (258, 1,048, and 4,054 actual usable points), and public T001
with x=`355fb5df7a905f16921eb`, y=`5900a390f42d290f1bbe`. The dimension-12
index must contain exactly 8,219,485 unordered pairs. The reference is
commit `2a0eb9b2e3fc6c24319754b2a4d77c2e01e1202a` with batch size 256.

Build release binaries from both source heads using the same compiler flags.
Run the small-base probe once per arm, then run the full-base probe in
baseline, candidate, candidate, baseline order. Preserve every exit code,
stdout, and stderr. Check curve/subgroup membership, pair counts, the exact
answer or exact miss, and group replay of any witness. Record build and query
nanoseconds separately, plus peak process RSS on the full base. Keep setup
and index build outside the target-dependent query time.

Retain the 2,048-element batch only if the full-base query median improves
by at least 10%, all small-base query medians regress by at most 10%, index
build medians regress by at most 10%, and peak RSS rises by at most 10%.
Otherwise revert the code and preserve this negative experiment. The host
is not isolated, so even a passing ratio is an exploratory stage diagnostic,
not a controlled CPU speedup or a complete F6/IC result. This four-summand
path has negligible natural relation coverage at dimension 12; no
one-target DLP or F4/F5/rho comparison is inferred from its timing.
