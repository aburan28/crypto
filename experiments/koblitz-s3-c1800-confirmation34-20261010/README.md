# Fresh 34-point N53 confirmation panel

Checked Sage generated **34 distinct N53 public target points**, T213–T246,
from the frozen SHA-256 scalar law. Every point has its own one-target workload
manifest and immutable point-file digest. [`PANEL_FREEZE.json`](PANEL_FREEZE.json)
fixes **136 fresh-process solves**: two C1200 reference and two C1800 combined
runs per target, with a seeded target shuffle and alternating ABBA/BAAB arm
order. Each process receives one public target after reusable preparation.

The planned threshold is a C1200/C1800 geometric online wall ratio of at
least **2.0** and a target-resampled two-sided 95% lower bound strictly above
**2.0**, with all 136 solves and all 34 independent scalar replays verified.
Execution status is **prepared: 0/136 process solves launched**. The sample size comes
from the previous N53 panel's empirical log-ratio variance; its 80% power
estimate depends on that variance carrying over to the confirmation host.

The [execution freeze](EXECUTION_FREEZE.json) binds the panel to the two
cross-built Linux x86-64 ELF executables (SHA-256
`3b5bbcab65333dc69c9ffb1acbaf41a3b5f4cd6256730641f590843eb0d8c3fb`
for C1200 and
`47223962d292361ed21bbf0ff577cfe3564c74a84773b16b3b467a5c1e24b8dd`
for C1800), their build receipt, and the runner, analyzer, and independent
Sage replay source hashes. [`run_panel.py`](run_panel.py) follows the frozen
136-run order in fresh processes, records raw failures and all five exclusive
online phases, and resumes only after a verified prefix. [`run_one.py`](run_one.py)
records each child's peak RSS and archives its full report, stdout, and stderr.
[`analyze_panel.py`](analyze_panel.py) retains incomplete status and computes
the target-resampled interval only after all 136 receipts verify.
[`replay_sage.py`](replay_sage.py) independently checks every scalar against
its public point. The offline preflight passed against both exact ELF files
without executing a target; the analyzer's empty-run probe returned unknown
timing and zero verified runs.

The comparison deliberately varies both actual factor-base size and target
extraction kernel: C1200 has 127,200 distinct subgroup-usable base points,
while C1800 has 190,800. The exact candidate IDs, source hashes, target
points, workload IDs, 14 preparation workers, one target worker, 600-second
per-process cap, and 32 GiB observed-RSS envelope are in the freeze. The
host isolation receipt remains pending. The runner currently expects a host
checkout containing the frozen point and workload files. The separately
prepared Modal T205 pilot transfers only its declared binaries and public
point bundle; the 34-target panel has not been uploaded or scheduled there.
A Modal guest can collect paired correctness and guest fidelity diagnostics;
the controlled-speedup gate additionally needs the
host-wide isolation and noise evidence specified in
[`ISOLATED_BENCHMARKS.md`](../../docs/ISOLATED_BENCHMARKS.md).

[`GENERATED_TARGETS.json`](GENERATED_TARGETS.json) retains the point and
workload file hashes. The generation script is
[`make_targets_sage.py`](make_targets_sage.py), and the checked Sage
[`runtime receipt`](sage_runtime_info.json) was saved before it ran. The
[`panel preparation script`](prepare_panel.py) verifies the point, workload,
candidate, source, prior-panel, and power-planning hashes before writing the
freeze. Generation and freeze preparation can be replayed locally with:

```sh
/Volumes/SSD990/cryptanalysis/sage --runtime-info > sage_runtime_info.json
/Volumes/SSD990/cryptanalysis/sage -python make_targets_sage.py
python3 prepare_panel.py
```

The commands above are intended from this experiment directory. On replay,
`prepare_panel.py` rejects an existing freeze; compare its archived hash or
write into a fresh directory.
