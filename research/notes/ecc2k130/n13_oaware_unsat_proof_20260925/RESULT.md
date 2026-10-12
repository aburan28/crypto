# Five fixed n13 O-aware PDP negative branches have externally checked DRAT proofs

The [pre-outcome protocol](PROTOCOL.md), exact #785 inputs, CaDiCaL binary,
vendored DRAT-trim source, five queries and caps were frozen in draft
[PR #790](https://github.com/aburan28/crypto/pull/790) at `a10007f`.
`FROZEN.json` SHA-256 is
`8941c739264bd0dc94eceb1d6d2efe3e108e13f6ee9f5ecec898f42456bd3cc4`.
The draft's hash-only CI and pre-outcome peer review passed before the first
selected proof run. The first and only five-case run began 2026-09-25
15:51:07 UTC and ended 15:51:14 UTC with `status=success`; every case and
both rejection controls passed. An independent fresh compilation of the
pinned checker then replayed all five committed certificates successfully.

| Fixed full-point branch | #785 assumption literal | CaDiCaL proof exit | External DRAT-trim result | Raw proof bytes | Solver wall s | Checker wall s |
|:--|--:|--:|:--|--:|--:|--:|
| Q0T0 | 676 | 20 | VERIFIED | 8,490 | 0.326 | 0.634 |
| Q4T0 | 629 | 20 | VERIFIED | 8,490 | 0.385 | 0.746 |
| Q4T1 | 504 | 20 | VERIFIED | 8,490 | 0.400 | 0.527 |
| Q4T2 | 985 | 20 | VERIFIED | 8,490 | 0.476 | 0.684 |
| Q4T3 | 690 | 20 | VERIFIED | 8,490 | 0.388 | 0.762 |

Each query byte stream reproduced its SHA-pinned #785 input: the exact #781
1,195,344-clause base plus one registered unit assumption. CaDiCaL 3.0.1
emitted text DRAT (`--no-binary`); the independent checker was compiled from
MIT-licensed upstream DRAT-trim commit
`2e3b2dc0ecf938addbd779d42877b6ed69d9a985`, source SHA-256
`d834b649f437e091597f5347f259b9f681087f89ca0844d0cee250a1a1a0c2ee`.
The fixed non-unit XOR control accepted its valid proof and rejected a mutant
with essential derived steps deleted. Q0T0's actual proof was also rejected
against the known-SAT Q0T3 query. The [retained first pre-outcome preflight
failure](PREOUTCOME_PREFLIGHT_FAILURE.md) explains why an earlier mutation of
only a redundant final step was an invalid checker control; it preceded the
freeze and every selected proof run.

These certificates establish UNSAT of **these five exact CNF branch queries**.
Q0T0 is a negative full-point torsion branch of an otherwise supported
projected Q0, so it does not make Q0 projected-negative. All four torsion
branches for Q4 are checked; together with #781's complete toy
encoding-to-point replay and #767's complete point oracle, Q4 is a
projected-negative **toy n13,m5** case. The DRAT-trim C implementation and
pinned CNF semantics remain explicit trust boundaries. The result does not
certify the other 22 #785 negative branches, direct S6/S7 encodings, n19 or
n131 PDPs, solver ranking, independent relation rank, a challenge logarithm,
full index-calculus cost or a rho crossover.

The five solver children used 1.976 seconds of summed wall, and the five
external checker children used 3.353 seconds; the wrong-input rejection check
used 0.950 seconds. These are separate stage costs, not an attack-time
comparison. All solver/checker children were below 60/180 seconds,
1 GiB RSS, 128 MiB raw-proof, 320 MiB temporary-workspace and 512 MiB
archive caps. The largest recorded cumulative child RSS upper bound was
161,300,480 bytes. Raw proofs total 42,450 bytes; deterministic gzip copies
total 11,672 bytes. The [first receipt](evidence/receipt.json), SHA-256
`5b051bc65788c8cd022b9e6440f778edd837ae3b8d3b3a67a8011dc5f4d9ae01`,
retains all raw proof hashes, commands, intervals, exits, streams and caps;
[EVIDENCE.md](EVIDENCE.md) gives the archive-only replay command.
