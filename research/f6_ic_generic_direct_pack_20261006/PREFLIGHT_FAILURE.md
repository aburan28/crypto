# First frozen run: declared-environment gate rejected the workload

The first freeze used implementation commit `8c635f185` and worker binary
SHA-256 `b2ae7023db05dc09bdd99f6948d7b959b5f5d29ad2590253138ed7b83f54fc88`.
All eight scheduled processes exited 2 before starting a target solve with
`undeclared algorithm environment override: KIC_F4_ACTIVE_MULTIPLIERS`.
No one-target timing or correctness result exists for those rows. Their
exact candidate/input manifests are in `preflight-freeze-v1/`; raw stdout,
stderr, exit status, timestamps, and file hashes are in
`preflight-failure-v1/` and `PREFLIGHT_FAILURE_SHA256SUMS`.

The worker intentionally rejects undeclared algorithm-changing environment
variables during exclusive measurement. The repair adds an explicit
`config.active_multipliers` Boolean and accepts exactly the matching
`KIC_F4_ACTIVE_MULTIPLIERS=1` setting, rejecting either a missing declared
override or an undeclared one. Both comparison inputs now declare the
setting. The changed worker, inputs, binary, and candidate identities are
frozen anew before any rerun.
