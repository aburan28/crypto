# N83 K3300 root11 canonical-selector screen

This bounded screen compares the ordinary 83-rotation canonical-root scan
with the existing `KIC_SPARSE_CANONICAL=1` selector. Both arms use the same
K3300/root11 factor-base construction, five-byte exact root slots, fixed
`2^30` table slots, 11-byte reconstructed-root states and Bloom disabled.
The sparse selector finds candidate rotations from cyclic zero gaps. Its
existing Rust tests compare key and shift against the scan exhaustively
through degree 13 and on sampled degree-83 words. The measured-arm verifier
checks the complete factor base and exact root-key count again.

The curve is `EC1N83Ckb1a1haacb69ee9f37`, polynomial-basis `GF(2^83)`
with modulus bitmask `9671406556917033397649557`, model
`y²+xy=x³+x²+1`, subgroup order `8569786107849059`, and cofactor
`1128547018`. Each arm must produce 547,800 usable points, 3,300 folded
columns, base hash
`a4ba33098d2776520b98a326585a2301f8ccd5bc005df7c1a70fbe8496c61d8b`,
903,852,075 root keys and final rank 3,300. Source SHA-256 is
`c9824acebda0360416c3c507594bb07c17ca798c32cbf8f84bba95ca2c44d5b1`.

Both fresh-process arms use the same 1,000 public points from
`research/n83_factor_base_refine_20261014/targets.jsonl`, SHA-256
`35e34151dcc1ad6048939560adf28fc7c1b1017833fe005a656c3f4e2bc07de8`.
The scalar sidecar SHA-256 is
`f8cbc90f7f592aae5014ca4f2693a3327a72e935970f22e8f02868e4aeed7f71`.
These points were used in an earlier base comparison, so this is a
secondary exploratory selector screen. A changed primary candidate needs
fresh same-point Tier 1 candidate/rho validation.

Run `canonical_scan` then `canonical_sparse`, with four rank threads, one
target thread, rank seed 29, a 1,200-second child cap per arm and a 14 GiB
sampled RSS guard. Compute `sparse_ms(i) - scan_ms(i)` for all 1,000
targets, partition them into 20 contiguous blocks of 50, resample the
blocks 10,000 times with seed `2026101508`, and use sorted 2.5% and 97.5%
endpoints. Advance the sparse selector only if both arms pass complete-base,
four-point-relation and scalar replay within the guard and the upper
endpoint is below zero. Preserve capped and failed statuses.

Record source, binary, fixture and host digests, complete setup and
process costs, peak RSS, state/table bytes, root-key count and every
target's online phases and probes. The online interval starts after
factor-base/index/log preparation and includes target query, point
decomposition, relation check, descent and recovery check.
