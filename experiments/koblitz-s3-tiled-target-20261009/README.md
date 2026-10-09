# Deterministic S3 target tiles

Fixed serial-order tiles make the N53 target query's S3 root-call count
repeat exactly on a fixed public point while 14 workers process disjoint
subranges. The eight-point T025–T032 size sweep verified all 64 one-target
executions and their recovered scalars; checked Sage independently replayed
all 64 and 462 distinct relation witnesses. Tile sizes 2,048, 4,096 and
8,192 gave exploratory immediate-parallel/tiled online ratios of 0.417×,
0.404× and 0.612×. The per-tile worker startup and the requirement to
complete each tile outweighed the saved root calls in this sweep, so none
advanced to the T001–T024 evaluation set.

## Pipeline and ordering

The exact curve is `EC1N53Ckb1h738a572b91f2`: polynomial-basis
`GF(2^53)` modulo `x^53+x^6+x^2+x+1`, `y²+xy=x³+1`, prime subgroup
order 21044858204113, 25,864 actual subgroup-usable factor-base points,
244 folded columns and 243 prepared relation rows. The
[parent freeze](../koblitz-s3-crosscurve-20261009/FREEZE.json) fixes the
public points, input law and workload IDs. The [tuning freeze](FREEZE.json)
fixes T025–T032, the exact implementation and each tile size; exact
manifests are under `identities/candidates/`.

| Arm | Candidate ID | Target tile size |
| --- | --- | ---: |
| Immediate parallel reference | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0he884cf865574` | rotating starts |
| Tiled 2,048 | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0h031b01b9610a` | 2,048 states |
| Tiled 4,096 | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0haba0a124f79e` | 4,096 states |
| Tiled 8,192 | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0h9f485af8e276` | 8,192 states |

For each tile, the algorithm partitions its consecutive serial offsets
among 14 workers. Each worker returns its earliest verified relation; the
smallest offset among workers wins. A tile with no relation is followed by
the next tile. This preserves the serial orientation law and first-witness
order independent of thread scheduling. All worker computations, startup,
tile barriers and draining are charged to that target's online interval,
which also includes target validation and scalar replay. Reusable base,
index and rank preparation occur before the online clock.
The [serial first-witness control](local/serial_control.json) ran the
group-lift source once on each of the eight points and found the same
relation indices in all 48 tiled executions.

## Matched size sweep

The [raw ledger](local/raw_tuning64.jsonl.gz) has one row per run keyed by
candidate, workload and run ID. Each of the eight public points was solved
twice by each arm in a frozen shuffled order. The [analysis](local/analysis.json)
checks all run identities and phase sums, retains failures, and checks that
the tiled root-call count is identical across repetitions. The
[Sage replay](local/sage_replay.json), with its checked
[runtime receipt](local/sage_runtime_info.json), verifies target relation
equations, prepared rows, rank and recovered scalars.

| Tile states | Exploratory immediate/tiled online ratio | Points above 1× | Immediate/tiled S3-call ratio | Tiled root-call A/A p95 abs log |
| ---: | ---: | ---: | ---: | ---: |
| 2,048 | 0.417× | 2/8 | 1.138× | 0 |
| 4,096 | 0.404× | 1/8 | 1.089× | 0 |
| 8,192 | 0.612× | 1/8 | 0.884× | 0 |

The local host's A/A wall-time ratios remain exploratory. The deterministic
root-call result isolates the cost of repeated worker creation and tile
completion as the next implementation target. A persistent worker pool can
keep the same first-witness rule with one spawn per target.

The tiled source SHA-256 is
`58dbc20c99f119c121db86dd8d79c90fab73dba858ad1132e6189a7e77dfdcd0`,
its local release binary SHA-256 is
`25864d4f6c14ada0b6a882a64459c2b66afb06260ec1d9c2717bded2ac9f0dc7`,
the plan SHA-256 is
`7256e134811196e896fd06429bf9ee9dcccad5399408815b230331d25b813428`,
and the uncompressed raw ledger SHA-256 is
`2b8a4a4f8e5f548f59fee0cf1d4d03d6548ed825da99521c3f18fd25edc0fafd`.
