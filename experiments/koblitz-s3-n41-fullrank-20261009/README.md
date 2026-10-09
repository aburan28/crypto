# Full-rank S3 direct recovery on N41

Completing the reusable factor-base log system lets a verified four-point S3
target relation recover its scalar from four signed Frobenius labels. The
N41 C400 configuration reaches rank 400 before target input and uses direct
recovery with zero target span checks in its T025 correctness control. It
retains polynomial roots and lifted points in the index and computes 41
inverse Frobenius target conjugates inside the online interval.

The [freeze](FREEZE.json) fixes the Koblitz curve
`EC1N41Ckb1hcbf4a37f6e66` over `GF(2^41)` with modulus `x^41+x^3+1`,
cofactor 4, and prime subgroup order 549756390943. The baseline has 244
signed-Frobenius columns, 20,008 usable points, a 16,777,216-slot root
table, and a rank-243 preparation threshold. The [C400 candidate](identities/candidates/IC1N41Ckb1fb32800PDP4rootRCsampleLAgaussTDdirectISO0h6023089a52e1.json)
has 400 columns, 32,800 usable points, a 33,554,432-slot table, and full
rank-400 preparation. Its retained polynomial coordinates add 209,920,000
bytes. Both use the same natural low-x construction and public curve; base
size and complete-log recovery are the declared algorithmic differences.
The exact bases are retained as [C244](local/n41_base_c244.json.gz) and
[C400](local/n41_base_c400.json.gz) encodings with digests in the freeze.

The frozen paired panel uses public N41 targets T025–T032, two repetitions
per arm per target, ABBA/BAAB order, relation seed 20260928, 14 reusable
relation workers, and one target worker. Each process loads one public point;
the online clock begins with target-dependent computation after reusable
preparation and ends after scalar recovery and verification. The five charged
phases are target query, PDP, relation check, descent, and recovery check.
All failures and timeouts remain raw ledger rows. The [runner](run_paired.py),
[analyzer](analyze_paired.py), and [checked-Sage replay](replay_sage.py) verify
the candidate identity, source/binary hashes, base digest, rank, public
point, recovered scalar, charged phases, and independent group relations.

The T025 C400 correctness control matched the public scalar, reached rank
400, used direct four-label recovery, and made 150 target S3 calls. Its
timing is excluded from the paired panel. A controlled CPU wall-time ratio
requires the repository's host-wide isolation receipt and noise gate; the
panel records exploratory timing and A/A repeats on this local host.
