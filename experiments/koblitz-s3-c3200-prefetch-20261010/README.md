# C3200 mixed quotient-root target prefetch

The C3200 target extractor now computes up to eight states' ordered S3 roots,
prefetches the corresponding mixed quotient-table home slots, then performs
the same exact lookup and group-lift sequence. The first verified relation
and all scalar checks remain in the target online interval. This changes the
target PDP memory-access schedule while retaining the C3200 factor base,
root-table representation, preparation rank, and C1200-first traversal.

The [frozen T025 control](CONTROL_FREEZE.json) binds source SHA-256
`56ed0e2167c9b5b1cc389b4d95ed7887f65bf17273d68911b9924ee4167008ac`,
release binary SHA-256
`54644a12e05ab3212dbad3dbf723f66dc913d1fb66bdaec135128570abb3dfae`,
the existing C3200 projective report, public point, runner, 14 preparation
workers, one target worker, 1,200-second cap, and 32 GiB observed-RSS cap.
Its pass criteria include the same first relation witness, 339,200 usable
base points, rank 3,200, exact mixed quotient table, the five exclusive
online phases, and recovered scalar. All 17 release example tests passed.

The control passed: rank 3,200, the same first relation witness
`[55035, 67241, 57152, 189185]`, independent Sage scalar replay, and
28,387,180,544-byte peak RSS. Its T025 target online interval was 1.450375 ms,
versus 3.036625 ms for the earlier projective control; target PDP was
1.353333 ms versus 2.939542 ms. The prefetch batch completed 512 S3 calls
versus 498 in the reference. These single-run timings are diagnostics on
an unisolated host.

The [frozen paired screen](SCREEN_FREEZE.json) compares this source with
the matched-rebuild C3200 projective reference on the eight previously used
T173–T180 public points, in alternating arm order. Both arms use the same
14 preparation workers, one target worker, 3,200-column rank, and 32 GiB
RSS cap. The primary screen statistic is the geometric mean of paired
single-target online ratios; every run must verify and fit the memory cap.
The screen advances only at 2×. Local CPU ratios remain exploratory until
host-level isolation and noise checks pass.
