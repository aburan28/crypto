# Frozen interface-admission protocol for the rotated m5/m6 solver corpus

Status: preflight only, frozen before a new solver process or benchmark. The
point-only Q targets, four torsion cosets, complete positive witnesses and
certified projected-negative targets are inherited byte-for-byte from merged
[#767](https://github.com/aburan28/crypto/pull/767). The #770 rational mask
reference, #774 affine S3 candidate audit and #777 finite-fibre theorem are
separate accepted evidence. This gate asks whether any *currently implemented*
S3-chain, direct-resultant or SAT-export path yields the same **complete
projected point PDP** on n13-m5 and n19-m6. It does not design a new exporter,
run a solver campaign, or estimate ECC2K-130 speed.

## Fixed inputs and admission rule

`INPUT.json` pins both curve arms, 16 Q labels through the #767 raw archive,
source/evidence SHA-256s, audited interfaces and a 30-second/256-MiB static
preflight cap. The canonical decision for each Q is an OR over its four
`Q+T` full-point targets, `T in E[4]`; a positive must yield a rational,
slot-valid signed point tuple adding to one exact full target, and a negative
means all four cosets are empty in the complete independent group-law census.
Do not treat a single empty torsion branch as projected UNSAT. The eight
positive and eight negative Q labels were frozen before any solver outcome.

An engine is admitted for timing only if an actual frozen input file exists
for **both** arities and all 16 Q labels, with one canonical variable/domain
map and explicit coverage of every required item in `INPUT.json`: distinct
rotated slot domains, rational factor lifts and both signs, all four torsion
translates, finite and O/inverse/identity branches, complete model-to-point
lifting, and sound SAT/UNSAT status. A direct S6/S7 algebraic zero is only a
candidate until rational signed point replay; a plain affine S3-chain zero
cannot prove full projected UNSAT without the exceptional branches. The
[#777](https://github.com/aburan28/crypto/pull/777) theorem proves affine
S3 finite-fibre soundness for sign-complete rational factor fibres, not
coverage of O prefixes or an exporter. The #774 n13 frozen target has the
exceptional-only x-mask `[0,0,0,2,1]`, even though that same Q also has
other affine witnesses; the control proves model-set incompleteness, not a
wrong target-level verdict on this particular 16-Q panel.

The preflight shall compare pinned source interfaces and prior reference
receipts, not infer absence from a repository-wide string search. Inspect the
current generic `build_decomposition_system`, binary Semaev functions, S5
SAT example, framework solver adapters and WDSat adapter. Enumerate local
binary availability, hashes and version output for CryptoMiniSat, MiniSat,
msolve, WDSat, Kissat and CaDiCaL on PATH, plus the prior #764 pinned
WDSat fixture path if present on this host. A missing local fixture on another
host is an inventory fact, not a source/corpus failure. Installed binaries alone are **not**
admitted arms: each needs a complete frozen encoding and a parser/verifier.
The existing S5 example has four factor points, and its fixture is not
silently retargeted to five or six distinct rotated slots. The source audit
is bounded to these identified interfaces; any new or overlooked complete
exporter needs its own explicit path and independent equivalence check.

If no interface meets every requirement, output `BLOCKED` with zero admitted
solver arms, the exact failed prerequisites, existing corpus SAT/UNSAT
certification and binary inventory. Preserve the raw receipt and stop: no
solver wall/RSS comparison is meaningful. If an eligible arm is found,
record it for a *new* preregistered timed PR; do not launch it under this
preflight. Do not turn the negative target certificates into a solver proof.

## Costs, replay and next gate

A cold Python child runs once with a 30-second process timeout and 256-MiB
peak-RSS acceptance cap. Record UTC interval, argv, exit status, stdout,
stderr, result SHA-256, binary path/version/hash, CPU/wall/RSS and any failure
file. The lightweight archive replay recomputes the pinned source/evidence
facts and deterministic admission verdict; it does not rerun #767's full
point oracle or any solver. The hash-only CI must pass before this preflight
result is read. Any source drift or failed check is preserved and not called
`BLOCKED` by a successful audit.

A follow-up exporter PR must freeze complete source equations/circuit bytes,
all exceptional branches and Q+T selectors, rational lift/domain constraints,
exact model-to-point decoding, variable order and engine-specific CNF/XOR/ANF
bytes before timings. Exhaustively compare its projected model set with the
#767 oracle on all 16 Q targets; include #774's exceptional mask as a model
control and an explicit no-affine/only-O-prefix planted control if needed.
Only then run paired cold SAT/UNSAT arms with pinned binary hashes, one-thread
flags, 15-second/2-GiB pilot caps from [#763](../ROTATED_M5_M6_SOLVER_ADMISSION_20260925.md),
full raw status/model/proof receipts, independent rational point witnesses,
and certified negative labels. FES needs a separately verified quadratic
encoding; n19's raw 88-bit S3 chain needs a multiword exporter. A solver-stage
result leaves full-DLP S and matched rho ratios unset.
