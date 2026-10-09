# Prelifted group-sum S3 index

The indexed S3 roots can be lifted to curve points during reusable setup.
Computing x coordinates of those points plus and minus a public target then
returns the same unordered S3 roots with one shared field inversion across
the two indexed points. This variant recovered all eight fresh public points
in a four-pair N41/N53 pilot; checked Sage independently replayed all **32**
one-target runs and **714** distinct relation witnesses. Its exploratory
incremental online ratios against the rank-243 prepared path were **1.14× on
N41** and **1.23× on N53**.

For curve points `P=(u,v)` and `Q=(t,w)` with `u != t`, let
`d = u + t`, `lambda_plus = (v+w)/d`, and
`lambda_minus = lambda_plus + t/d`. For the frozen `a2=0`, `b=1` curves,
the two S3 roots are
`lambda_plus^2 + lambda_plus + d` and
`lambda_minus^2 + lambda_minus + d`.
The [Rust candidate](../../examples/s3_partial_prepared_grouplift_targetspan.rs)
prepares one y coordinate for each indexed root using `z^2+z=x+1/x^2`,
`y=x*z`; it retains the scalar S3 path for exceptional or nonliftable states.
Every one of the **2,440,976 N41** and **3,155,408 N53** states in the frozen
subgroup index had a valid lift, so the measured target runs used the group
path throughout. A unit test checks the unordered root sets on both curves.

The [freeze generator](freeze_grouplift.py) records source SHA-256
`433eb4368b42bd51894d0ea1f41681914a9190f90a0703c01c4a566069aa5ed2`,
local binary SHA-256
`82285b41518bb1f9fe5bc8f6c94807f131a4510fdd7252aeec82526e91bcd91b`,
and [`FREEZE.json`](FREEZE.json) SHA-256
`ee815bcecc399f6a655e5585e024e586cd71c4213ffb22e0898b8887457da18c`.
The [candidate manifests](identities/candidates/) specify the exact curves,
actual base sizes, folded columns, rank 243, source hashes, and group-lift
policy. The [matched-curve plan](../koblitz-s3-crosscurve-20261009/README.md)
supplies the same 200 frozen one-target workload IDs; this pilot used indices
1–4, two runs per arm, ABBA/BAAB order, 14 workers, and a 180-second per-run
cap.

| Curve | Prepared candidate ID | Group-lift candidate ID | Actual base points |
| --- | --- | --- | ---: |
| N41 | `IC1N41Ckb1fb20008PDP4rootRCsampleLAgaussTDdirectISO0h0f38cf1ef60a` | `IC1N41Ckb1fb20008PDP4rootRCsampleLAgaussTDdirectISO0hd7645060ef00` | 20,008 |
| N53 | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0hbd278c9cf111` | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0hcc3dbaef9275` | 25,864 |

| Curve | Verified runs | Prepared / group-lift online ratio, geometric mean (minimum) | Target PDP ratio | Cold ratio | Online root-query ratio | Prepared / group-lift A/A p95 absolute log ratio |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| N41 | 16/16 | 1.136× (0.673×) | 1.264× | 0.493× | 1.000× | 1.221 / 0.295 |
| N53 | 16/16 | 1.232× (0.970×) | 1.231× | 0.972× | 1.000× | 0.836 / 0.859 |

The [paired analyzer](analyze_pilot.py) retains every run and checks target,
candidate, workload, binary, scalar, preparation rank, online phase sum,
valid index lifts, and relation witnesses. The local wall ratios are
exploratory: all A/A p95 values exceed the 10% gate
`abs(log(repeat1/repeat2)) <= log(1.10) = 0.0953`, and no target achieved a
2× incremental online ratio. Group-lifted state preparation also raises the
N41 cold time on this host. A host-wide isolation receipt is required for a
controlled CPU wall-time comparison.

The [raw pilot archive](local/raw_pilot4.tar.gz), SHA-256
`6cdaa3e93ed0337f13ca3981dcfd97f6177e85b341e1cf12a1d8e05fc9aa2e60`,
contains the 32 JSON reports, stdout, stderr, ledger, and analyses. The
[analysis](local/analysis_pilot4.json) records each point's timings and root
queries. The [direct Sage replay](local/sage_replay_pilot4.json), SHA-256
`400f1d7b3a6b36d91e1c1f1bd059cb32a160e67781c0f559b2f2686e21adeedf`,
uses the checked [runtime receipt](local/sage_runtime_info.json). Its
[equation-equivalence receipt](local/sage_equivalence_pilot4.json) also checks
the new path against the earlier Sage-replayed rank-243 relations.

Reproduce with the repository's checked Sage launcher:

```sh
CARGO_TARGET_DIR=/private/tmp/s3-prepared-target-20261009 cargo test --release --offline --example s3_partial_prepared_grouplift_targetspan
CARGO_TARGET_DIR=/private/tmp/s3-prepared-target-20261009 cargo build --release --offline --example s3_partial_prepared_grouplift_targetspan
python3 experiments/koblitz-s3-grouplift-20261009/freeze_grouplift.py
python3 experiments/koblitz-s3-grouplift-20261009/run_paired.py --first 1 --count 4 --binaries /private/tmp/s3-prepared-target-20261009/release/examples --output /private/tmp/s3-grouplift-paired4-replay
python3 experiments/koblitz-s3-grouplift-20261009/analyze_pilot.py --ledger /private/tmp/s3-grouplift-paired4-replay/raw.jsonl --output /private/tmp/s3-grouplift-paired4-replay/analysis.json
env LD=/usr/bin/ld /Volumes/SSD990/cryptanalysis/sage --runtime-info > /private/tmp/s3-grouplift-paired4-replay/sage_runtime_info.json
env LD=/usr/bin/ld /Volumes/SSD990/cryptanalysis/sage -python experiments/koblitz-s3-grouplift-20261009/replay_sage.py --ledger /private/tmp/s3-grouplift-paired4-replay/raw.jsonl --runtime-info /private/tmp/s3-grouplift-paired4-replay/sage_runtime_info.json --receipt /private/tmp/s3-grouplift-paired4-replay/sage_replay.json
```
