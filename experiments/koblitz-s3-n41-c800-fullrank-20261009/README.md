# N41 C800 full-rank S3 holdout

The 800-column low-x factor base solves its complete log system during
reusable preparation. A verified four-point S3 relation then recovers one
public target scalar from signed Frobenius labels, with the 41 inverse
target conjugates charged inside the online interval. The T025 correctness
control reached rank 800, verified the scalar, used direct recovery, and
made 26 target S3 calls. T025 was used to select this configuration and is
excluded from the timed holdout.

The [freeze](FREEZE.json) fixes Koblitz curve
`EC1N41Ckb1hcbf4a37f6e66` over `GF(2^41)` with modulus `x^41+x^3+1`,
cofactor 4, subgroup order 549756390943, and 24 new public targets
T033–T056. The [C800 candidate](identities/candidates/IC1N41Ckb1fb65600PDP4rootRCsampleLAgaussTDdirectISO0hb8875aea0063.json)
has 65,600 subgroup-usable base points, 800 signed-Frobenius columns, a
134,217,728-slot S3 table, and 839,680,000 retained indexed-coordinate
bytes. The baseline has 20,008 usable points, 244 columns, a 16,777,216-slot
table, and a rank-243 preparation threshold. Both bases use the same
natural low-x construction, with exact [C244](local/n41_base_c244.json.gz)
and [C800](local/n41_base_c800.json.gz) encodings retained for Sage replay.

The [runner](run_paired.py) fixes two repetitions per arm on each point in
ABBA/BAAB order, relation seed 20260928, 14 target-independent relation
workers, one target worker, and a 120-second cap per run. Each process
receives one public point. The online interval starts with target-dependent
work after reusable preparation and ends after scalar recovery and check;
its exclusive phases are target query, PDP, relation check, descent, and
recovery check. The [analyzer](analyze_paired.py) retains failures and A/A
repeat spread, and the [checked-Sage replay](replay_sage.py) verifies the
encoded bases, group relations, matrix rank, and recovered scalars.

The target labels were frozen after the C400 pilot and before C800 timing.
The factor base size and full-rank recovery are declared differences from
the baseline. Local wall-time ratios remain exploratory unless the CPU
host-wide isolation and noise gates are met.
