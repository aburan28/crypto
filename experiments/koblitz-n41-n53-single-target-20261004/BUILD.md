# Frozen producer build and local environment

The implementation commit is `6773fb1b11483ecdf0a05fa61d031f213e08d048`
on top of crypto main `fd15d65586b5e1966678dd726a8c308ecca2ecf8`.
The [hash manifest](BUILD_SHA256SUMS) pins the three release example binaries,
their source files, and the copied [Cargo lockfile](Cargo.lock.frozen). The
build command was:

```sh
cargo build --offline --locked --release \
  --example koblitz_orbit_dlp_fast_online \
  --example koblitz_rho_fixture \
  --example koblitz_one_target_replay
```

The local host is a MacBook Pro Mac16,7 with an Apple M4 Pro (10 performance,
4 efficiency cores), 48 GB RAM, macOS 26.6 (Darwin 25.6.0, arm64), and
Homebrew `rustc 1.93.1` (commit `01f6ddf7588f42ae2d7eb0a2f21d44e8e96674cf`,
LLVM 21.1.8). The host does **not** provide an auditable exclusive CPU
partition, fixed frequency, IRQ routing or a host-level noise receipt.
Its wall times can only be exploratory. `RAYON_NUM_THREADS=1` and
`OMP_NUM_THREADS=1` are set for each run. The time wrapper archives process
wall/CPU separately from the producer's charged in-process intervals; the
producer records peak RSS with `getrusage`. An initial n41 rho solve under
`/usr/bin/time -l` succeeded but the wrapper itself exited one after sandboxed
`sysctl kern.clockrate` failed. Its raw output is retained under
`runs/n41/preflight_time_l_failure/`, and its exact Q is retained as
`targets/n41.jsonl`. Later n41 runs reuse that Q without substitution. The
runner now uses `/usr/bin/time` without `-l`.

Before target generation, the new producer and verifier were built with
`--locked`, the archived n37 result replayed successfully, and the new cold
phase path was exercised on that archived n37 point. The n37 cold-control
receipt reported `cold_phase_verified: true`, 518 actual base points, seven
folded columns, rank seven, and independent scalar replay. The n37 control
is a functionality check, not a new n41/n53 outcome or an isolation claim.
