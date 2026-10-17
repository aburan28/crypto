# Direct inverse-Vélu fiber oracle, frozen protocol V1

Frozen before implementing or running the new oracle. This is a public,
single-target ECDLP comparison on the selected different-`j` cell
`p=1237`, `F_(p²)=F_p[s]/(s²-ν)`, `r=1201`, cofactor `h=1275`, source
`y²=x³-21x+37`, kernel `T=(3,1)` of order three, and the existing
degree-three Vélu quotient. The source/destination `j` values are
`1066/727`. The source subgroup generator, 64 distinct public targets,
destination images, old full-image base, and independent checker answers
are fixed by `../order1000/BATCH_1237_1201_V2.public.json` and its
separate `.checker.json`. The public input SHA-256 is
`24f2c56a9ba56f91c954f38b7752edf24fee5ed027b49b1ac0747771945099c4`.
The solver must not read the checker input.

## Hypothesis and mathematical check

For a destination affine point with abscissa `u`, a source preimage
abscissa `X` satisfies the monic inverse-Vélu polynomial

`G_u(X)=X³-(u+6)X²+(6u+21)X-(9u+32)`.

This follows from the explicit map's x coordinate
`u=2(y²+1)/(X-3)²-X-12` and the source equation, for `X≠3`.
The experiment tests whether a direct finite-field polynomial root
oracle can identify the first eight projected image-base columns
without enumerating the entire source group. A rational root is only a
candidate: lift it on the source, evaluate the full Vélu map, and accept
only an exact destination point match. Reconstruct all three rational
preimages by adding `0,+T,-T`; choose the same lexicographically first
source witness as the frozen source census. Handle the destination
identity and exceptional source points explicitly. Finite checks do not
replace a proof of the general fiber statement; a Lean proof map must
cover the polynomial identity and mark any curve-level bridge still open.

## Frozen execution

Scan destination points in ascending `ExtCurve::key` order. Factor the
degree-three polynomial over `F_(p²)` using `gcd(G_u,X^(p²)-X)` and a
deterministic quadratic-character split with shifts `c=0..31`. A split
cap failure, repeated root, bad lift, map mismatch, or RSS/child cap is
recorded and stops promotion; no fallback source census may be used by
the oracle arm. Count destination points scanned, polynomial calls,
field operations, split attempts, rejected points, exact map checks,
group operations, CPU and peak RSS. Existing field counters omit
`sqrt_mod` internals, so any `S_group` remains a lower bound until those
native steps are priced.

Run three deterministic rounds and prefixes `1,4,16,64` with five arms:
`inverse-velu`, `full-image`, `source-bsgs`, `rho-negation`, and
`rho-negation-aa`. This is 60 cold child processes and
`3*(1+4+16+64)*5=1275` independently checked scalar outputs if all
finish. Use the existing `8r` relation cap, 32 rho restarts of `20r`
steps each, `ceil(sqrt(r))` BSGS baby steps, 300 seconds per child,
and 512 MiB peak RSS. Preserve the exact public target order and seed
formula from `different_j_trend_batch.rs`. The old full-image base is
an equality control, not input to the inverse oracle; compare raw and
projected columns, canonical source witnesses, relation support keys,
rank, relation witnesses, scalar outputs, and final source and
destination scalar replays. Do not pass known answers to any solver.

The generic signed-rho floor in counted group operations is
`sqrt(pi/(2*2))` in `S_group`; the measured same-target reference is
`rho-negation`, with the 64-target old order-1201 value
`S_group=4.071` as a preregistered context value, not a new control.
The complete comparison table has one `S_group` column for every arm
and separate `F_p` multiplication/inversion, polynomial, CPU, and
memory axes. `S_group` includes every counted setup, map, relation,
solve, and replay group operation divided by `N sqrt(r)`. It is a
lower bound on portable `S` until the omitted native categories and
host conversion are priced under `../PRICE_PROTOCOL_V2.md`.
Descriptive CPU ratios are reported only with a same-host rho A/A
control; an isolated wall comparison still requires the separate L2
gate. Keep cold `N=1` and reused setup `N=4/16/64` separate.

Success for this bounded run requires exact base reconstruction,
complete verified relation recovery and scalar replay on all 64 public
targets in each round, a native independent checker for all child
outputs and relation witnesses, matched controls, a source-linked
report/graph/PDF, and all raw failures retained. A lower group-count
ratio alone is a structural diagnostic, not a speed claim. If the
oracle changes the base or stalls within the frozen cap, preserve that
result and inspect the mathematical exception before changing the
protocol. A second-host exact replay and at least four spaced orders
remain separate follow-on gates for a scaling claim.
