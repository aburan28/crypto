# Continue only the registered SAT jobs after an exporter-manifest key error

The first dispatch of the [frozen panel](panel.json) produced all four
explicit-point wide-S4 exports at exit code 0 under the registered exporter
watchdog. The original runner then marked **every** control `INVALID_EXPORT`
with `KeyError: 'mitm'`: its verifier read `manifest['mitm']`, while the
pinned exporter writes the documented key `direct_meet_in_the_middle`.
Every raw manifest has `native_sat.status: not_run_in_export_process` and
`direct_meet_in_the_middle.status: not_run_in_export_process`; the exporter
did not run either internal solver. The original runner returned normally
after retaining all four failed rows, and **no CryptoMiniSat process was
started**. This is an adapter-schema failure, not evidence about SAT or the
existence of point decompositions.

The original 42-file [stage-A archive](stage-a-evidence.tar.gz) is committed
unchanged; SHA-256
`2b1fa7abb093c94e741cebecdce35ea8fd59dbc0fb301e9aac9f1b95e0b195ef`.
Its `summary.json` hash is
`b75c9e2c2f621bdef499e2f441d772a4435d93d1b1bfc49278ad8adafb68778c`,
and its preregistered runner hash is
`fc555329cbe4fb299d03fac68fecf22d5a6c29abd9743bf5fefff25ddf015f79`.
The exporter executable hash is
`4da7f5781da1aba2f76d6b9ce8f6e4ea43845f2465cf8edf0612615931d0629a`.
All original exports, manifest bytes, process metrics, stdout/stderr and
build receipts remain available inside the archive. The original summary
continues to describe what the first runner decided.

The versioned [continuation runner](../../continue_cms_s4_controls.py) must
first verify the full stage-A archive and every original row, source/build
receipt, binary hash, point and exact label. It adapts only the incorrect
manifest-key lookup by checking the real
`direct_meet_in_the_middle.status` and then passing that value to the frozen
validator under its expected local key. It rechecks all declared source
exports before any SAT call. If any other mismatch occurs, it stops without
starting a solver. Then it executes the **four previously registered**
CryptoMiniSat jobs once, in the same trial order, from the retained XOR-DIMACS
files. It uses the original pinned binary, parser, process meter, thread,
random seed, conflict cap, model count and 120-second watchdog. It does not
rerun the exporter, select new queries, increase limits or overwrite the
original output. Raw SAT streams, model validation and independent group
lift are retained in a new output directory.

This repair is committed and pushed to the same PR **before** continuing
the four solver jobs. The claim boundary remains binary-bound correctness
diagnostics only. A first-model source SAT assignment without a group point
witness remains nonlifting; it cannot qualify a complete source-bound IC
solver. Any timeout or failure remains an outcome and is not retried.
