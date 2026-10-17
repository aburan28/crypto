# Frozen secondary local batch inputs

The 128 public points in `batch_targets.jsonl` and verification-only scalars
in `batch_answers.jsonl` were generated and hashed before the fine-grid
batch. They are disjoint from the 24-point prospective holdout and earlier
N83 fixtures. `batch.py` and `verify_batch.py` are byte-for-byte snapshots
of the scripts frozen in `local_batch_plan.json` and
`local_verify_plan.json` in the parent directory. Their absolute paths
identify the original local artifact and binary, where full raw outputs are
retained; run the scripts from that artifact rather than this source copy.

The batch starts only after the K2752, K2880 and K2944 single-target
arithmetic checks in `../single_target/` have passed. Its five arms are
K2752, K2816, K2880, K2944 and K3072. The preregistered local selector
requires independent replay of all 128 answers, bases and four-point
relations for every arm, then chooses the lowest geometric mean online
interval with a lower-K tie break. This batch selects an arm for the separate
24-point Tier 1 confirmation; its shared preparation and ordinary host
measurements are secondary diagnostics.
