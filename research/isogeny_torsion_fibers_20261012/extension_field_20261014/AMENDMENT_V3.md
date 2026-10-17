# Paired-run scope amendment before DLP execution

The [frozen fixtures](FIXTURES.jsonl) put the prime-order subgroup in the
anti-Frobenius eigenspace, so its points have `x` in the base field. The raw
line `x=s*t` is outside that subgroup and has a large cofactor component.
A direct raw pair-sum relation must accidentally cancel that component;
the raw-pair arm records this exact control. Projection by `[h]` gives a
checked subgroup relation because `[h](P_i+P_j)=[h]R`, and its work is charged
from base construction through recovery.

The [paired spec](RUN_SPEC.json) adds a complete cofactor-projected pair-table
solve, a `chi` orbit **index** on those exact projected sums, and the explicit
degree-three isogeny control. The `chi` key is only a candidate index: every
returned relation is checked as an exact source-subgroup point equality.
This is a bounded, fully timed extension-field DLP comparison of those
implemented arms. It does not implement the `chi` polynomial-system oracle or
the direct Semaev resultant specified in [the original protocol](PROTOCOL.md).
Those original oracle and isolation requirements remain open and are not
weakened by this added discriminating run.

The one-target cold interval includes line lifting, cofactor projection,
pair-table construction, every query, relation recovery, linear algebra,
and verification. Same-target rho includes its walk construction and scalar
check; the six-class rho also charges eigenvalue discovery in the measured
whole-process CPU interval. The primary same-unit table uses thread CPU work
per calibrated source addition divided by `sqrt(r)`, because the native
coordinate transforms, allocations, and root extraction are all executed
inside that interval. This is an L0 host-dependent unit, not a portable
primitive count. A second table retains exact counted operations and
unpriced-native caveats.
