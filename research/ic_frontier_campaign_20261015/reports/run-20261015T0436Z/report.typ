#set page(paper: "a4", margin: 2cm)
#set text(size: 9.5pt)
= IC frontier campaign — Phase 0 pilot run report

*Run.* 2026-10-15, sessions started 04:36–04:52 UTC. Nine `ecbench` sessions; specs `pilot/spec-*.json`;
binary SHA-256 in `pilot/ecbench.sha256`. Host: repository cloud container, Linux x86_64, 4 vCPU, 15 GiB;
`--cpus none --allow-busy`, so every record is *L0*: operation counts only, no wall time admitted.
Protocol: `PROTOCOL.md` (draft, not frozen).

*Workloads.* Public hash-to-subgroup targets, train seed 20261015100; one or two targets per arm.
Curves: Koblitz $E_a$ at $n = 13, 19, 23$; prime: `prime_search` 20-bit, PKM ladder 24/28-bit,
`find_a3_curve` 32/36/40-bit.

*Requested and completed.* 65 executions; 57 verified, 8 exhausted: negfold with full-rank linear algebra at
$n = 13$ (2) and both prime IC arms with a fixed 64-point base at 32, 36 and 40 bits (6). No timeouts, errors or crashes. Every record is in `pilot/TABLE.md`.

#figure(image("pilot_s.svg", width: 100%),
  caption: [Charged $S$ per arm and curve, as recorded (single targets; open marker = lower bound,
  × = exhausted). Sessions ECBS1h8da761a6e558, ECBS1h4e3c36fa57bd, ECBS1h7d1741682d1e, ECBS1h8a85c79c72d0,
  ECBS1hf28feb0fe98c, ECBS1ha6179c6e7e0d, ECBS1hfda2c8befe2c, ECBS1h94b78342b569, ECBS1hda3b226c750b.])

*Observations (no comparison is made).* Every rho, BSGS and kangaroo record is a lower bound. IC negfold and
F4 records are complete at $n = 13$ and $n = 23$ (pinned calibrations) and lower bounds at $n = 19$.
`ic.shared_rank` lists three unpriced categories without counts. *Correction:* an earlier version of this
report said `ic-compact` and the prime IC arms do not emit their unpriced counts. They do, in each phase's
`native` map; the error was in the summary script (`pilot/FINDINGS.md`, finding 2). One inherited-F4
execution at $n = 23$ took about 6 minutes (charged $S$ = 85,071). Interpretation and protocol
consequences: `pilot/FINDINGS.md`.

*Raw records.* `pilot/sessions/<name>/records.jsonl`; replay audits not run for the pilot.
