# Frozen third-cell selector for the distinct N=64 gate

Status: preregistered before the p=401 construction or pairing run.
This cell extends the earlier batch protocol without changing either
historical input or result. Candidate curve: `E/F_401: y^2=x^3+1`, with
`F_(401^2)=F_401[s]/(s^2-nu)` and `nu` the first positive quadratic
nonresidue. Check `#E(F_401)=402` by exhaustive x/Legendre counting, then
derive and check `#E(F_(401^2))=161604=402^2` from the trace recurrence.
Use `r=67`, test primality, verify `67^2` divides the extension order and
the source generator has exact order 67. If any check fails, stop and retain
the failed screen.

Select the first finite F_401 point in ascending x and then ascending y
whose `[6]P` is nonidentity; set `G=[6]P`. Set
`k=1+SHA256("pairing-batch/2026-10-15/p=401/r=67")[0..8]_BE mod 66`
and `Q=[k]G`. The solver input contains only the public coordinates and
group parameters; the scalar fixture is separate. Form all 66 distinct
nonidentity targets with the original Q first and the remaining scalars
shuffled by `StdRng(seed=0x4241544348000000 XOR 401)`. Freeze the files,
hashes, H-search seed, and native source before timed runs.

Use the existing reduced Tate pairing, one independently selected H in
`E(F_(401^2))[67]`, one alpha, a reused multiplicative BSGS baby table,
and one beta/finite-field DLP/scalar check per target. Certify the degree-two
field, H's order, alpha's exact order 67, and pairing bilinearity on 1, 2,
and 3 multiples before comparison. This source subgroup lies over F_401;
the order-three curve automorphism is not assumed to preserve it. Price
same-target charged negation rho, with the V3 step cap and 32-retry schedule;
record the absence of a six-class subgroup action explicitly.

Run three rounds and cold cumulative prefixes `N=1,4,16,64`, plus `N=66`
as a complete-group endpoint. Charge field construction, every H candidate,
Miller and final exponentiation, the finite-field table and solve, failed
attempts, scalar replay, memory, thread CPU, wall, and rho restarts. Keep
the original limits: 4096 candidate x, 512 acceptable H, 10^7 field
multiplications, 512 MiB RSS, and 15 min thread CPU per row. Independent
native replay must verify every scalar against the same public Q. A failure
or cap exhaustion remains a censored outcome. Compare cold prefix totals
and per-target marginal work separately; N=64 is an amortized result.

The reference floor for one generic DLP is `sqrt(pi/4)` counted operations
per `sqrt(67)`. Raw CPU is descriptive L0 shared-host evidence until every
field and group primitive is converted to the same counted unit.
