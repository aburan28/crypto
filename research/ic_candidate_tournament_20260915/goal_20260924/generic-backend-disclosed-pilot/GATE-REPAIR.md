# Recover the predeclared stage-B gate from raw incomplete reports

Schema v3 produced the two n17a1 stage-A worker reports once and did not
execute stage B. The original runner's [summary](RESULT-v3-stage-a.json)
marked both jobs as process failures because their exit code was 2. The pinned
CLI deliberately returns code 2 for a valid `status: incomplete` report when
`max_trials=1` has not solved the entire DLP. It still writes the full query,
PDP, base, matrix and phase record. The raw stdout SHA-256 values are:

| Arm | Raw stdout SHA-256 | PDP outcome | `unsupported` | Queries |
| --- | --- | --- | --- | ---: |
| F4 | `2f24cd7914634665c03dbf53102276bbcc4a6790c491c193381b3a129509a550` | witness | false | 1 |
| F5 | `bde8532abe74a1c073600413349db16d4c318fdefed181a1d3bf5b65dd34dc46` | witness | false | 1 |

Both jobs finished within the 60-second cap, with sampled RSS far below the
8 GiB threshold. Independent build binding, exact standard-subspace base
reconstruction, natural query-law/group replay and complete generic stage
audit all pass on each retained report. These results satisfy the *original*
stage-A condition of one audited natural PDP attempt and
`stats.unsupported=false`; they do **not** satisfy complete-DLP qualification.
The original runner summary SHA-256 is
`ff58ad86bebd13000347c5971f2c29da42374cd4ef1700cee4035be28e3e4a23`.
It remains unchanged as an account of that runner's decision.
The controlled build record SHA-256 is
`fe6378068da024ffcfdd4d5d23d7d2c1221284ffd6fbe810ca28510177188a93`.
The committed [raw stage-A evidence](RESULT-v3-stage-a.json) includes the
complete two worker reports, their jobs and process outcomes, the original
summary, build record and source manifest. The panel and raw jobs use
`algorithm_seed=2026092918`; the frozen protocol prose says `2026092917`
in one sentence by mistake. No stage-A input or measurement was changed to
resolve that typo.

The continuation script is versioned before further execution. It rechecks
the registered schema-v3 panel, the controlled build and binary, both raw
stage-A jobs and their independent audits, and that no stage-B job ran. Only
then does it run the eight **already preregistered** n19/n23/n31 jobs once,
with the same source, solver configs, disclosed points, seed and limits.
It accepts a valid incomplete worker report at exit code 2 for a stage
diagnostic while retaining null complete-solve costs. It never reruns the two
stage-A jobs or overwrites their original summary.

Run `continue_generic_backend_disclosed_pilot.py` from the tournament
directory with `--source-root` pointing at the pinned source checkout,
`--original` pointing at the retained stage-A bundle, and `--out` at a new
directory. The script refuses a changed panel, original summary, build,
stage-A stdout, stage-B directory or source checkout before it creates any
output. It writes a corrected stage-A audit receipt separately from the
original result and journals each stage-B job. The pass rule remains all ten
preregistered cells dispatched and independently audited; any failure stays
in the final table.
