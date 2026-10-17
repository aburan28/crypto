# N83 cooperative target-worker cancellation screen

The N83 four-summand indexed decomposition already supports multiple target
workers, but its current policy joins every worker after each has found a
relation. The new opt-in `KIC_TARGET_CANCEL=1` policy stops other workers at
their next state boundary after the first worker finds a relation. It keeps
the exact root index, factor base, rank system, and independent group check.
The returned witness is the first completed valid relation; worker scheduling
may change which valid relation is returned.

Run three source-identical K992/root11 policies sequentially on the same 128
public N83 points, with fresh processes and empty caches:

| Policy | Target workers | Cooperative cancellation |
| --- | ---: | --- |
| `serial` | 1 | off |
| `parallel_join` | 2 | off |
| `parallel_cancel` | 2 | on |

Use the exact field, curve, generator, and subgroup in the existing N83
canonical record. The public fixture is the previously sealed
`research/n83_root11_finegrid_holdout_20261017/local_batch/batch_targets.jsonl`
(SHA-256 `72486b7d708023838db1104a4a5e75a04b915d6f6d4013057c934fb31000aa87`),
with its answer sidecar used only by the independent checker. The compiled
source, binary, runner, verifier, commands, environment and raw output each
receive SHA-256 receipts before or after their corresponding execution.
Each policy has a 1,200-second process cap. Preserve timeout, nonzero exit,
missing output, and verifier failures as censored outcomes.

For every completed policy, replay the full base, all 128 four-point relations
and target scalars, the five target phase costs, and the group check. Require
identical exact base hashes across policies. Compare per-target online wall
intervals and full cold-process intervals separately. The 128-point run is
an ordinary-host shared-setup stage diagnostic. A fresh-point Tier 1
same-point IC/rho comparison follows only if the cancellation policy yields a
verified improvement on this screen.
