# Preregistered charged n19 phase-selector successor

Status: **design freeze only; no solver outcome is authorized by this note.**
The exact [phase census in PR #807](https://github.com/aburan28/crypto/pull/807)
passed its predeclared cap-8 support screen: the spread order reached 128,751
of 130,873 subgroup targets with 276,497 **free archive-oracle probes**. Its
four-phase prefix (121,419; 253,355 probes) did not dominate the fixed
four-beta comparator (121,851; 249,827). Five betas reached 126,265 in
258,849 probes. None of those probe counts prices an implicit PDP. This
successor asks whether the three-column phase policy still helps after
complete six-point solver queries, certified misses, base-log acquisition,
and target recovery are charged. The exact literals and source hashes are
in [INPUT.json](INPUT.json); `check_protocol.py` checks the pins and case
split without opening phase outcomes.

## The one admissible PDP representation

Build a **new complete n19, m=6, d=2 exporter**. The selected factor in
each of six *distinct* rotated slots is exactly one of that slot's seven
rational **finite** full points. Every slot includes the 4-torsion point
`(0,1)`, which is not O; no factor choice is O. Use one-hot selector groups, without
cross-slot ordering or a shared `F0` domain. Map selectors to canonical
`(o,x,y)` flags/coordinates and chain five full-point additions to the
exact full target `tau^s(Q)+T`, where `T` runs in the frozen four-torsion
order. Every prefix can be O; copy, inverse, self-inverse x=0, doubling,
and distinct-x addition must all be represented. There is no affine-only
escape, existential `x` root without a rational sign lift, or shortcut from
`[4]` equality to full-point equality.

The preferred implementation reuses the full-point XOR/AND relation of
[#802](https://github.com/aburan28/crypto/pull/802) and ordinary Tseitin
DIMACS gate of draft [#804](https://github.com/aburan28/crypto/pull/804),
but those are respectively a local relation and a **single-edge** exporter.
Neither is already an S6 input. Their merge/release gates must complete,
then the six-slot chaining, one-hot links, target pins, model decoder,
proof/checker runner, and independent verifier must be committed in a
separate implementation PR. If a different complete encoding is used,
freeze its exact equations, variable order, semantics and byte hashes in
that PR before any solver outcome; the arm and target policy below remain
unchanged. A partial exporter does not enter timing.

An n19 five-edge circuit has 11 distinct point roles (six factors, four
prefixes, one final target) and five slope words. At 19 bits, those are
`11*(2*19+1)+5*19 = 524` primary point/slope bits before 42 one-hot
selector bits, Tseitin auxiliaries or constants. The old 88-bit affine S3
chain already exceeded its `u64` mask. This exporter must use arbitrary-width
IDs and reject `DIMACS id > 2^31-1`; a reported 524- or 566-bit layout is
**not** a claim of a compact CNF. No `m!` reduction is sound across the
different factor slots.

Use one pinned CaDiCaL 3.0.1 binary as the primary cold solver; its current
known SHA-256 is frozen in `INPUT.json`. Use text DRAT from that exact
binary and the separately compiled upstream `drat-trim` source at commit
`2e3b2dc` to certify every UNSAT branch. A SAT result needs exit 10, exact
status, a complete consistent Boolean model satisfying the independently
parsed CNF, exact selector/point/auxiliary decoding, and an independent
full-point sum equal to `tau^s(Q)+T`. UNSAT needs exit 20, a nonempty proof
under the proof cap, and a successful external checker on the exact query
bytes. An unknown, timeout, OOM, missing or invalid model/proof is
`CENSORED`, never UNSAT. Preserve every raw attempt, including failures.

### Semantic admission before the 64-target run

1. Freeze the implementation source, exact input and binary SHA-256 values,
   host/toolchain, resource limits, commands, and archive parser in its
   implementation PR. Its hash-only CI and independent peer review must
   pass at the **exact head** before any solver call. The source phase PR
   and #802/#804 must be merged; verify the actual merged hashes and re-freeze
   if they differ from today's draft. Do not silently substitute a binary.
2. Independently enumerate all `7^6=117,649` labelled beta-3 factor tuples
   at n19 and the four other frozen beta arm factor sets. For each tuple,
   derive its four-torsion target by a separate bit-serial/Fermat group law,
   construct a complete satisfying circuit/model assignment, and validate
   its CNF. Compare the resulting beta-3 projected support with #767's
   complete 62,389-point archive and each other beta histogram with #769.
   Independently evaluate the circuit relation on every branch type,
   including O and inverse prefixes; a matching hit count alone is
   insufficient. The n13 five-slot control must include #774's exceptional
   mask `[0,0,0,2,1]` for `Q+O=(7256,3272)` and all 32 #785 branch
   labels, with complete model lifting and external proofs for all negative
   branches. Then run the eight frozen n19 Q from #767 (four planted, four
   exact negatives), each with all four torsion branches. Every positive
   needs a lifted point witness; every negative branch needs external proof.
   Do not use the archive to choose solver target order or skip a query.
3. An exact mismatch, incomplete branch, unproved negative, width/cap breach,
   or checker failure stops admission. Publish the raw failure and a corrected
   **new freeze** before retry. A local #804 one-edge pass alone cannot
   release the six-slot campaign.

### Fixed policies and same-Q work

The point-only target file from #779 fixes 64 `ho-000` through `ho-063`
coordinates without scalars. The disjoint 256 `tr-000` through `tr-255`
training points are processed in that literal file order. Training labels
may form the right side of a verified relation but cannot be supplied to
the PDP solver or target recovery child. All arms use the same stream, fixed
torsion order `O,(0,1),(1,0),(1,1)`, one fresh solver process per branch,
and the same caps. Run arms independently from cold state; do not borrow
base logs, query results, clauses or proof caches across arms.

The primary phase arm uses **only beta 3** and spread phases
`[0,5,10,15,1,6,11,16]`, in that order, for both training and holdouts.
Its first four phases are a predeclared prefix diagnostic. A separate cold
beta-3 control uses phase 0 alone. The comparator arms use beta order
`[3,338435,303097,464276]` and the same list plus `42605`, each at
phase 0. A beta or phase is called negative only after all
four torsion branches have checked UNSAT proofs; stop at the first verified
SAT witness, and do not query later branches for that attempt. If an earlier
branch is unknown, retain it and continue the fixed schedule; a later SAT
can still yield a valid point log, but a final miss remains unresolved unless
every attempted branch was proved negative. Never use #807 first-hit bytes
to choose which phase or beta to call.

For training, process each Q through the arm's fixed schedule until the
first verified witness, form the signed `[4]F0` row, and update exact rank
modulo `q=130873`. Stop only at full rank (expected column ceilings 3,
12, 15; actual distinct columns must be measured) or after all 256 Q.
The phase row from target `tau^s(Q)` must be multiplied by
`lambda^(-s) mod q`, with `lambda=41811`, to represent `[4]Q` in the **same
three** signed beta-3 columns. Check each factor point, torsion shift,
cofactor projection, row group identity, training equation, column log
point, and exact matrix rank independently. Rank deficiency blocks that
arm's IC recovery; it is not repaired with sealed holdout scalars.

On each holdout, follow the same fixed arm schedule until a verified SAT
row or cap exhaustion. Recover `k = 4^(-1) * row dot column_logs mod q`,
verify `[k]H=Q` by independent group arithmetic, and only then open the
sealed label. For a fully checked negative policy, report all proof cost;
for any unresolved miss, report censorship and all attempted work. A
predeclared automorphism-aware `±tau` Pollard-rho reference must solve the
**same 64 Q** on the same host and be checked by `[k]H=Q`. Every IC arm
falls back to that matched rho implementation for any Q it has not solved,
so the final comparison has 64 verified logs per arm and charges its
failed/negative PDP work **plus** the fallback. Keep the primary stage
coverage, first-witness, complete-negative, and fallback counts separate.

Use one primary rho seed per case:
`little_endian_u64(SHA256(domain || "/rho/" || case_id || "/0")[0:8])`.
Two more seeds `/1` and `/2` are frozen secondary repeat controls. The
implementation freeze must pin the rho code and prove its quotient-walk
automorphism semantics, cycle handling and scalar recovery against the
ordinary group law. It must not inspect the sealed holdout scalar. Where
possible interleave matched arm runs by case and repeat; report the exact
schedule. A single run permits a **cost diagnostic**, not a runtime speed
claim. A claim of runtime improvement needs paired cold reruns and a 95%
paired interval excluding parity. A cryptanalytic `S` or rho crossover
claim also needs a defensible measured common *operation* unit for SAT and
rho; CPU/wall alone are reported as platform diagnostics with `S=null`.

## Caps, cost ledger, and decision

The frozen limits are in `INPUT.json`: 2,000,000 DAG nodes, 256 MiB per
query CNF, 1 GiB proof, 60 s/2 GiB export, 120 s/4 GiB solver, 300 s/4
GiB proof checker, 60 s/1 GiB independent model replay, and 72 h/16 GiB
for one full 64-case arm. Process-tree RSS and external monotonic wall
enforce caps. Count file I/O and proof checking in the charged wall/CPU
ledger. Stop an arm on its whole-arm cap and retain partial evidence; do
not normalize censored calls to zero. The semantic admission panel has a
separate 8 h/16 GiB cap. No cap is relaxed after seeing a status.

The cost table has beta-3, spread-phase K8, four-beta, five-beta and rho
rows; K4 is a labelled **prefix diagnostic** of the K8 arm with shared
training, not an independently cold attack row. For each arm report input/build, one-time base and
log training, every branch export, SAT, checked UNSAT, unknown/censored
attempt, model lift/point verification, rank solve, holdout descent,
fallback rho, and total CPU/wall/RSS. Put exact point-recovery count,
complete negative count, and 64/64 final correctness beside cost. Count
probes only as a separate structural column. For 64 targets publish both
**cold total** and warm `total/64`; never present a free-oracle mean as
solver cost. Derive operation-equivalent `S` and ratio to matched rho only
if a common-unit calibration is separately validated and frozen.

The phase path remains a priority for the next size only if it completes
the semantic and proof gates, attains full rank, recovers every solver-hit
Q correctly, and its charged 64-Q total including fallback is lower than
**both** four- and five-beta totals on equal 64/64 verified workloads. A
rho crossover requires a validated common operation cost with
`speedup = C_rho / C_IC > 1`, `S_IC = C_IC / sqrt(q)`, and
`S_rho = C_rho / sqrt(q)`, with a paired runtime interval if wall speed is
also claimed. If phase has no charged advantage, or caps censor any
negative needed for a decisive result, report a negative or censored n19
decision; retain a measured component bottleneck as a possible engineering
follow-up, but do not promote the phase selector from oracle probes alone.
Nothing at n19 extrapolates by itself to ECC2K-130: a later n37/n41/n53
and n131-relevant iteration must re-freeze its own matched costs and scaling.

The implementation PR must archive every exact CNF/proof/model hash, raw
stdout/stderr (including empty files), solver/checker exit, command, binary
SHA, UTC interval, wall/CPU/tree RSS, independent replay status, selected
full-point witness, row, rank and verified scalar. Its post-outcome CI must
reconstruct source hashes and factor domains, parse every archived CNF,
recheck every SAT model against clauses and group law, externally recheck
every UNSAT proof, recompute all rows/ranks/logs and the 64 final point
equalities. It may cite #807 as an independent support oracle only after
replay, never as a substitute for a negative proof or measured cost.
