# Exact retained K=600/K=2,048 mixed-orbit pair

The previously stored primary `public_x_hash`, seed `2026100801`, K=600
signed-Frobenius object is an **exact ordered prefix** of the retained K=2,048
envelope. Its 99,600 decompressed point rows match the first 99,600 envelope
rows byte for byte. The curve, subgroup, policy, seed and closure match;
accepted candidate indices, representatives and pre-projection source points
also match as ordered prefixes. The remaining envelope suffix has 1,448
signed orbits and 240,368 point records. This adds a second verified
small/envelope pair beside K=1,182/K=2,048. No relation search was run for it.

The K=600 object was read from its existing
`s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/a0/objects/b4dda45e6c3edfca673db6f3f82c16cce8c860f5e62c092155258df399a0363e.jsonl.gz`
key. The [S3 head response](s3-head.json) reported 3,025,798 bytes. The
downloaded compressed bytes have SHA-256
`3e85c9e419fcddd3c9faed3c84470bf96fcad9178e48e1d73629369335841ecf`;
their BLAKE3 digest was independently recomputed as
`b4dda45e6c3edfca673db6f3f82c16cce8c860f5e62c092155258df399a0363e`,
agreeing with the frozen [panel manifest](../../pilot-01/manifest.json) and
[upload/download receipt](../../pilot-01/upload-receipt.json). The envelope
download matches its prior audited SHA-256
`aacb62784dee0c9f5e06328dbb85a4b6a1c7243e778328da82900e7284f3a74a`.
The [machine screen](screen.json) binds both manifest SHA-256 values, the
original panel upload receipt, the S3 head response and the prior 15,000-second
search audit. [screen.py](screen.py) replays all headers and rows from the
exact local compressed objects without using a solver or changing S3.

For five summands, the exact unordered point-multiset count is
`C(B+4-j,5-j) × C(L+j-1,j)` with B=99,600 small and L=240,368 residual
points. The count divided by the nonidentity subgroup size is 0.033785426
for 5+0, 0.407661103 for 4+1 and 1.967593178 for 3+2. Their
uniform-target Markov hit-probability ceilings are respectively those first
two ratios and **one** for 3+2. The 3+2 ceiling is vacuous; it neither
establishes a relation nor predicts solver success on the frozen public
target. Six-summand counts and exact fractions are also in the machine
screen. The [small exact controls](test_screen.py) check disjoint-domain
multiset arithmetic and the threshold crossing.

Compared with the already admitted K=1,182/K=2,048 pair, K=600 cuts the
small-base column count from 1,182 to 600 but enlarges the residual orbit
pool from 866 to 1,448. Those changes affect matrix size, SAT domains and
large-prime graph closure in different directions. There are no natural
partial-yield, graph-cycle, rank, sparse-linear-algebra, individual-log or
matched-rho measurements for this pair, so it cannot be selected as the
minimum-total-runtime factor base. The earlier conservative local charge
remains 15,000 seconds; this read-only object comparison made no new retained
N83 SAT call.
