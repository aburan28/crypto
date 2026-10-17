# Timing correction before qualified comparison

The first complete executions used the committed solver at
`8906100efa`. `RUN_V1.jsonl` contains all 84 run rows and four calibration
rows, but `/usr/bin/time -l` exited 1 after the solver completed because
`sysctl kern.clockrate` was denied in the sandbox; `RUN_V1.stderr` retains
that infrastructure failure. The direct native rerun `RUN_V2.jsonl` exited
0, and `REPLAY_V2.json` independently checked 75 recovered scalars and
retained nine exhausted raw-pair rows.

The V1/V2 driver started its CPU and wall timer before serialising each
`RunReport` to JSON. The IC report is much larger than a rho record, so
those whole-row timings do not price equal method intervals. Their exact
group-operation, relation, rank, answer and replay data remain useful.
Do not use `RUN_V1` or `RUN_V2` CPU/wall ratios as a method comparison.

The V3 driver stops both timers immediately after the solver returns and
before `json!(report)` or any output. It retains the same fixture hash,
three rounds, seeds, arms, caps, factor bases and algorithms from
`RUN_SPEC.json`. The six-class rho interval still includes eigenvalue
discovery. Calibration and output serialisation are outside every arm.
`RUN_V3.jsonl` and its separate replay will be the timing basis after
successful execution; V1/V2 are not overwritten.
