# Host and run conditions

- Apple M4 Pro, 14 hardware threads, 48 GB RAM, macOS (Darwin 25.6), `cargo 1.93.1`, release profile.
- Run 1 (`bench.md`, `bench.json`): 2026-10-07/08 local time, command
  `prime_ecdlp_fast_bench --baseline-ic-bits 20 --baseline-s4-bits 16 --rho-max-bits 56 --ic-max-bits 32 --seeds 3`.
  The host was shared with two other sessions compiling and running experiments:
  load average 87–110 on 14 cores at the start of the run, 50–67 at the end,
  ≈ 20 GB of swap in use. Treat sub-second single-thread rows as noisy (see
  the (†) rows in `RESEARCH_PRIME_FAST_ECDLP.md`); multi-second rows are
  conservative.
- Run 2 (`bench2.md`, `bench2.json`, if present): sections `micro,rho` only,
  `--rho-max-bits 48 --seeds 3`, taken afterwards at lower load.
- All instances are public synthetic known-answer curves: the committed a=−3
  ladder rungs (16–28 bits) and `find_a3_curve` generations (32–56 bits,
  parameters printed inline in `bench.md`). Every recovered log was checked
  against the planted scalar before being reported.
