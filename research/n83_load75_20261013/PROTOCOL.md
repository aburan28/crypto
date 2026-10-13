# N83 packed-root 75% load: frozen comparison

## Question and selection rule

The opt-in `KIC_PACKED_ROOT_TABLE_LOAD75=1` requests
`ceil(2*regular_states/3)` entries from the existing twofold, power-of-two
allocator. This can keep K2202 in 8 GiB of packed slots under the source-loop
bound `regular_states <= 83*K^2`. The conditional numeric capacity proof is
`../n83_dense_root_20261013/Load75Capacity.lean`; the source-to-state bound
is recorded as a separate proof obligation there.

Compare three complete prefix-base pipelines on the same N83 public points:
K1798 with the existing dense request, then K2000 and K2202 with the 75%
request. Select the arm with the lowest mean target probes over all 64
verified panel targets, among arms completing rank and all targets within
900 seconds and 8 GiB of root-table slots. Include every timed-out/failed
arm in the results. Keep online wall times exploratory on this shared host.

## Frozen inputs and execution

- Curve: Koblitz a=b=1 over F_2^83, polynomial low terms [0,2,4,7],
  subgroup order 8569786107849059, cofactor 1128547018, generator
  (6614244651398184860292015,8868260410029249391539529).
- Factor-base construction: deterministic prefix `construct:83:1:K`,
  cofactor projected, signed-Frobenius closure. K denotes effective columns;
  actual distinct usable point count and base hash are measured.
- Public point file `targets.jsonl` is a byte copy of the frozen
  N83 dense-root holdout with SHA-256
  `e9ed7ddb0aeafdbfc7354d67279a8712572eea9d2aee3896f99633bdc2e937cc`.
  The first line is the primary unseen point in `primary_target.jsonl`;
  scalar answers remain separate in `answers.jsonl`.
- Rank seed 29, 12 rank threads, one target thread; each process capped at
  900 seconds. Primary order: K1798 dense, K2000 load75, K2202 load75.
  The 64-point secondary panel runs in reverse order. Each panel arm dumps
  its full point-defined base for independent replay.
- The single-target same-point rho control is the already frozen
  `fold-L64-DP6-E0` receipt in `../n83_dense_root_20261013/runs/rho_primary`.
  Its 60-second cap, seed 2026101301, point and correctness certificate
  are unchanged. The new producer never receives scalar answers.
- Capture process status, wall/RSS, exact output, base count/hash, state/key
  counts, table capacity/bytes, probe count, exclusive target phases, rank,
  and independently replayed scalar and four-point relation certificates.

## Code and binary identity

- Parent frozen dense-root result commit: `21fda0911dbeacbf2a0cc42246f012c6bb2ff78e`.
- Producer source SHA-256:
  `d4bf0f34d9d57eb26ed795a610be4baa469b3835455d69ff902ef46a20036fdd`.
- Producer release binary SHA-256:
  `fc4b62da8e5dc7980c3e2a27e774166d594b5661641833252e8ec8ae6894fa40`.
- Built with `CARGO_TARGET_DIR=/Volumes/SSD990/llm/tmp/n83-dense-root-cargo-target cargo build --release --locked --offline --example koblitz_orbit_dlp_fast`.
- Twelve focused release example tests passed. The 75% table lookup
  behavior has a direct 768-entry/1024-slot first-witness control.
- Local Apple M4 Pro, 48 GiB RAM, macOS 26.6, shared CPU. Host-wide
  isolation is unverified, so controlled speedup remains unknown.

`run.sh` checks the source, binary and public-point hashes before every
arm. Each run directory is write-once. Raw statuses remain immutable. The
first outer timer failed after its producer run; [AMENDMENT.md](AMENDMENT.md)
records the preserved failed attempt, corrected `primary_v2` timer, and
the peak-RSS wrapper frozen before the secondary panel.
