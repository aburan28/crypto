# Run report: prime root-finding pair-table size

Session `ECBS1h382073f2a6ad` executed the committed
[`specs/prime20-b32-vs-b64.json`](../../specs/prime20-b32-vs-b64.json)
from source commit `2273c71ee1d3adddf47d4e6a2d75f06669d43739`.
The exact command was:

```sh
RAYON_NUM_THREADS=1 target/release/ecbench run --wait --allow-busy --cpus none --spec research/ecbench_rootfinding_tuning_20261013/specs/prime20-b32-vs-b64.json --out research/ecbench_rootfinding_tuning_20261013/sessions/prime20-b32-vs-b64
```

The host capsule records macOS arm64, 14 logical CPUs and environment class
`ECBENV2hd82681268e96`; the binary SHA-256 is
`be956553a5a2cd1227e8f1add8dd24e11dffba5068c799234b55e44ffc3a1b61`.
The runner earned L0 on every execution. It wrote 128/128 planned records,
all verified; 96 were measured after warmup, 24 per arm. There were zero
failures and zero timeouts. The full replay receipt at [`../../replay-all.json`](../../replay-all.json)
reports 96 identical measured replays and no problems.

The recorded mean cold GAE values are 22,094.250 for B32, 18,343.208 for
B64, and 971.917 for matched negation rho. The saved paired comparison gives
B64/B32 = 0.830225 (95% target-cluster interval 0.697026–0.961546), 24
pairs across eight public targets. These are observed counted costs; native
prime work is unpriced and wall values are descriptive L0 only. The chart
source is [`../../render_chart.mjs`](../../render_chart.mjs), its values are
[`../../summary.json`](../../summary.json), and the graphic is
[`../../cost.svg`](../../cost.svg). [PDF report](report.pdf).
