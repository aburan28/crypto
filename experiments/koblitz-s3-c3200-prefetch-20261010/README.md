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

The control has not yet been run. If its correctness and memory gates pass,
the next measurement pairs this source with the C3200 projective reference
on frozen, previously used development points before drawing on new targets.
Local CPU ratios remain exploratory until host-level isolation and noise
checks pass.
