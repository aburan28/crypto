# Fixed-base n53 rank-restart pilot

This experiment holds the n53 compact-orbit curve, 23,320-point K220 factor
base, four-summand S3-root decomposition, pivot-guided rank policy, final
linear algebra, and development public point fixed. It compares unbounded
rank searches with searches restarted after 200,000 or 400,000 support probes.
The cap applies to rank collection only. Every abandoned attempt contributes
its query, decomposition, relation-check time, and probe count.

`PROTOCOL.md` is the preregistration copied from cryptanalysis PR #600.
`FROZEN.json` binds the exact protocol, producer source and executable, build
lockfile, base, target point, three candidates, three workloads, and resource
envelope. Candidate manifests are `candidate_*.json`; `workload_*.json` fixes
one rank seed each. The actual factor-base point count is 23,320 and the
folded matrix has 220 columns. The base digest is
`7af2460c8b5a2c29f9d1aa7fefecbcde3a6ce761293d0dab0cecfc3bc980b973`.

The producer was built on macOS arm64 (Apple T6041, Darwin 25.6.0) with
Homebrew Rust/Cargo 1.93.1 and LLVM 21.1.8:

```sh
CARGO_TARGET_DIR=/private/tmp/crypto-compact-target cargo build --offline --release --example koblitz_orbit_dlp_fast_online
python3 research/notes/ecc2k130/n53_fixed_base_rank_restarts_20261010/freeze.py --binary /private/tmp/crypto-compact-target/release/examples/koblitz_orbit_dlp_fast_online
python3 research/notes/ecc2k130/n53_fixed_base_rank_restarts_20261010/run_pilot.py --binary /private/tmp/crypto-compact-target/release/examples/koblitz_orbit_dlp_fast_online --dry-run
```

The recorded source SHA-256 is
`4dd3c5e7823b453ad1fa9046d9be68a04bddc79343ea29576dd1f32b0bc66299`;
the executable SHA-256 is
`3c03013e4a12d5eac74e1737c5549558ddf734aa9a324426914778765da0f51e`.
The executable is a local build artifact; the source and lockfile are tracked.
`run_pilot.py` checks those digests before execution, rotates the nine cells,
caps each process at 60 seconds and 16 GiB observed peak RSS, and preserves
the first result of every cell. It writes each process's raw files and
independent rank/target replay receipts under `runs/pilot/`.
`PILOT_RESULT.md` and `PILOT_ANALYSIS.json` give the paired audit and the
preregistered 400,000-probe cap selection for the held-out step.
`HELDOUT_GENERATION.json` pins the new public-point rule and paired seeds.
Its first execution stopped at a contradictory CLI backend assertion before
generating Q; `inputs/generation_attempt1/` retains the empty fixture, stderr,
and failure receipt. `HELDOUT_GENERATION_AMENDMENT.json` pins the minimal
strong-backend gate repair and unchanged hash seed before the second attempt.

The final cap=1 n13 control in `controls/n13_cap1_seed12/` reached rank 2/2 in five
attempts, three of which were capped. Its trace accounts for all five probes;
independent rank and target arithmetic replay passed. `verify_control.py`
checks its source and executable identity, rank trace, and exclusive phase
sums; `final_control_receipt.json` records PASS. The published scalar 7 in
that control is a correctness fixture. The earlier formatting-only build and
its control remain under `controls/preformat_identity/` and
`controls/n13_cap1_seed12_preformat/` with their original source and binary
digests, separate from the final producer.
Pilot input is the published point
`[7960849849661793,7443722527872608]` without a scalar supplied to the
producer.

Pilot wall times are exploratory on this shared host. The protocol selects a
cap only among candidates completing full rank and independently verified
target recovery on all three rank seeds, by the lowest median complete cold
cost. The six-seed held-out comparison and same-point rho pairing have their
own preregistered gate in `PROTOCOL.md`.
