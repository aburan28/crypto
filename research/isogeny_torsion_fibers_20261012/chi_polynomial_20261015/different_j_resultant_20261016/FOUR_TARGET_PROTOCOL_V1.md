# Four held-out targets and complete relation solves

Frozen after exhaustive χ exactness and before target generation or
relation-run implementation. This implements Gate 2 of
[`FOUR_GATE_PROTOCOL_V1.md`](../../portable_order_trend_20261015/inverse_velu/FOUR_GATE_PROTOCOL_V1.md).

Use the committed order-229 different-`j` public fixture with SHA-256
`f6c3126b4d95e3c457c9e552471bc2cdff322de0e0a9de9fbba701d7c05897ba`.
Enumerate all nonidentity `[k]G`, `1<=k<229`, and sort by the existing
`different-j-order-batch-target-v1:13:<point-tuple>` SHA-256 domain,
then by point tuple. Verify the first 64 public source/image points
against `PUBLIC_V1.json`; freeze entries 65–68 (zero-index positions
64–67) as the four new source/image targets. The solver-public file
contains no scalar or answer-bearing field. A separate checker file
contains only the four scalar correspondences and the public-file hash.
Hash both files and version the source/binary freeze before any solve.

Use the eight full-image source witnesses in the parent fixture,
cofactor-project them on the source and carry the eight signed columns
through the degree-three quotient. This gives every relation arm the
same eight columns under the explicit map. The χ arm uses the checked
chord resultant and roots from the V1 exactness gate; S3 source and
S3 quotient arms use direct third-summation-polynomial root solving;
the inverse-Vélu arm constructs its image base directly; the full-image
arm uses the source census. Each arm must check that its selected
projected columns match the common frozen base before a paired-yield
comparison. A mismatch is retained, with paired-yield unset.

Run χ, source S3, quotient S3, inverse-Vélu, full image, source BSGS,
signed source rho and rho A/A on each of the same four targets in
three rotated rounds. Relation arms use seed
`(0x20261015 XOR 211) + round*0x9e3779b9 + target_position*0x01234567`
with wrapping arithmetic, `8r=1832` trials per target, rank nine for
complete solve, and the same signed columns. Rho has 32 restarts of
`20r=4580` steps. All child processes have 300-second wall and
512 MiB peak-RSS ceilings. Charge setup, every polynomial and failed
root, map evaluation, query, rejected candidate, row operation,
recovery, scalar replay, CPU, wall and memory. Reuse across four
targets is secondary to each cold first target.

An independent Rust checker must replay every recovered scalar on
source and quotient, every relation witness, the common base, the
random streams and all phase/counter totals. Do not rank incomplete
or censored children as zero-yield. A complete portable `S` needs the
separate primitive calibration gate; process CPU and wall on a shared
host are descriptive. Identical pair support is the mathematical
control; the discriminant is verified useful rank per complete work.
