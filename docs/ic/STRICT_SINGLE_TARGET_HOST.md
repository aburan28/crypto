# Strict isolated one-target Koblitz IC/rho run

The paired adapter in `tools/ic_single_target_isolated.py` freezes one fresh
public point, runs compact-orbit IC and strong rho on that point, independently
replays both recovered scalars, and emits the internal online interval and
retired-instruction count to the strict isolated benchmark service. The
service records the physical-host preflight, paired order, raw failures, noise
checks, code hashes, and full stdout for each arm.

## Frozen instance and local correctness rehearsal

The manifest's default point law uses SHA-256 domain
`kic-n61-isolated-one-target-20261009-v1`, role `isolated`, and seed `0`,
followed by the protocol's half-trace lift and cofactor projection. The fixed
`K=400` comes from the earlier disjoint tune; the new point did not select it.
The input to both producers is the public point alone.

| Field | Value |
| --- | --- |
| Curve | Koblitz `a=0`, field degree 61, subgroup order `162888033982417` |
| Public `Q` | `[1382975434415487392, 1140017945195175241]` |
| Point hash in service manifest | `3f1b3745a8f9a4e305b95a1e08c71f90ea94e0b51e132f7fca2b4d9adf75fc6b` |
| IC factor base | 48,800 subgroup points, 400 signed-Frobenius columns |
| Independent scalar replay, both arms | `47476952596908` |
| IC online interval | first target query through scalar replay; 45.066542 ms on the local rehearsal |
| Rho online interval | first target-dependent walk operation through scalar replay; 233.141209 ms on the local rehearsal |
| Local resource and counter status | macOS ARM64, one producer thread, ordinary shared host; hardware instruction counters unavailable |

The local rehearsal checked the four IC relation points against `Q`; its five
exclusive IC target phases summed to 45.066542 ms. The rho record used an
empty distinguished-point table and its walk/collision and replay phases
summed to 233.141209 ms. These local times are diagnostics. The paired
physical-host result is produced by the procedure below. The adapter reports
`independent_replay_ms` separately from each producer's online interval.

## Physical-host procedure

Prepare a physical Linux host with the isolated cgroup v2 partition, separate
housekeeping CPUs, fixed frequency, memory policy, IRQ routing, and root-level
audit described in `cryptanalysis/docs/ISOLATED_BENCHMARKS.md`. Permit the
producer's `perf_event_open` user-instruction counter. Keep the two repositories
at pinned commits, and compile both Rust examples from the same `crypto`
checkout. The following paths are the host layout used in this example:

```sh
cd /workspace/crypto
cargo build --release --locked \
  --example koblitz_orbit_dlp_fast_online \
  --example koblitz_rho_batch_ks_strong_online

python3 tools/ic_single_target_isolated.py manifest \
  --ic-binary /workspace/crypto/target/release/examples/koblitz_orbit_dlp_fast_online \
  --rho-binary /workspace/crypto/target/release/examples/koblitz_rho_batch_ks_strong_online \
  --cgroup /sys/fs/cgroup/benchmark-isolated \
  --cpus 4-5 --execution-cpu 4 --mem-nodes 0 \
  --output /workspace/isolated-bench/n61-one-target.json

python3 /workspace/cryptanalysis/scripts/isolated_bench.py \
  probe /workspace/isolated-bench/n61-one-target.json
python3 /workspace/cryptanalysis/scripts/isolated_bench.py \
  --queue-root /workspace/isolated-bench \
  submit /workspace/isolated-bench/n61-one-target.json
```

The CPU and NUMA IDs above must match the host's verified topology. Start the
strict service on housekeeping CPUs before submitting. Its `probe` result must
pass, and its completed job must retain `preflight.json`, `manifest.json`,
`runs.jsonl`, `pairs.jsonl`, `summary.json`, and both raw stdout files. The
manifest requires a positive online retired-instruction count for each arm;
the service rejects a failed replay, a mismatched target or scalar, a missing
counter, an incomplete pair, changed artifacts, or a failed host/noise check.

For the primary result, copy the service receipt into the candidate run,
construct its exact candidate/workload/run identities from the retained base
and source hashes, and run the `vs_rho` claim checker. Preserve all raw failure
rows and the IC preparation costs beside the paired online interval. Promote
`rho_online_ms / IC_online_ms` only after the strict receipt and claim check
both pass.
