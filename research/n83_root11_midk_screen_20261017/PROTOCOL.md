# N83 root11 factor-base size refinement near K3072

This screen measures the target-online cost of the original five-byte root
index for K3088, K3104, K3114 and K3136. The previously proved allocation
boundary keeps K3114 within 2³⁰ root slots; K3115 crosses into 2³¹ slots.
K3136 samples the larger table. The original layout's proved key bound
permits K through 3596. K3072 is the frozen lower control and K3584 is the
independently replayed upper control.

| K | Root slots | Root table bytes |
| ---: | ---: | ---: |
| 3072 control | 2³⁰ | 5 GiB |
| 3088, 3104, 3114 | 2³⁰ | 5 GiB |
| 3136, 3584 control | 2³¹ | 10 GiB |

Use the sealed 128-point shared-setup fixture at
`research/n83_root11_finegrid_holdout_20261017/local_batch/batch_targets.jsonl`
and its answer sidecar. The IC source SHA-256 is
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`;
the already built release binary SHA-256 is
`faaac22aa6228c8d486875d6947fa654417d1122191fbe0eb9bde59bb45eb98f`.
`plan.json` freezes the exact paths, code and input hashes and root11 flags.
Each K starts a fresh process with empty caches and a 1,800-second cap.
Run in ascending K and stop after the first producer or independent replay
that fails. Do not retry or discard a failed raw arm.

For each completed arm, independently reconstruct all 166K factor-base
points, check full rank K, replay all 128 four-point group relations and
scalars, and reconcile the five exclusive target phases to the target
interval. Report geometric target-online time, median, mean probes, full
cold process, and root state/table memory. Pair each target interval with
the frozen K3072 interval and report the geometric ratio and descriptive
95% interval. K3584 is a separate upper control and does not get rerun.

This is an ordinary local-host secondary shared-setup discovery screen.
The descriptive minimum across the six measured K values is a candidate
for a newly sealed, fresh-point Tier 1 same-point IC/rho confirmation;
it is not an online factor-base selection from this screen alone. Preserve
all censored arms and the known K3584 orchestrator failure receipt in their
original archives.
