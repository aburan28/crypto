# ISO-1 theorem and larger-prime continuation (2026-10-08)

The contribution for this round is the generalized norm-one conductor
theorem and the exact absolute-Frobenius quotient. The original requested
range remains every prime from 37 through about 200, with every trace and
a necessary-and-sufficient class criterion. Completed diagnostic censuses
through 47 remain evidence; p = 53 is now complete and validated; p = 59 is running, followed by 61.

| Requirement | Execution and evidence |
| --- | --- |
| Publish the finished report | Existing PR #1588; extend the same PR. |
| State a theorem | `THEOREM.tex`, compiled PDF, explicit halving and quadratic-order proofs. |
| Continue to larger numbers | Absolute-orbit census queue for primes 53 through 199; preserve one status and receipt per prime. |
| Every trace | One Hasse candidate CSV row per trace `t = 2 (mod 4)`, including zeros and nonordinary rows. |
| Necessary condition | Prove it at every size; generalized PARI controls test the formulas. |
| Sufficient condition | Keep open until an existence proof and certified census labels support it. |

## Frozen enumeration and controls

Use the field and group arithmetic extracted from commit
`62becf9572fe74cbe8b3d8cebee3bf8a240708a1`; the source provenance and
the focused Cargo package are under `census_runtime/`. The two-point
Hasse-interval counter remains randomized. The package compiles the
retained arithmetic kernel independently of the unrelated library errors.

The `--derive-twists --absolute-orbit-quotient` mode enumerates the
normalized alpha values, selects the least encoded parameter among
`lambda^(+/-p^j)` for `0 <= j < 6`, and weights its trace by the number
of distinct orbit elements. For `p > 3`, assert
`point_counts = (p^4+3p^2+8)/12` and
`weighted_representatives = 2p^4+2p^2`.

Before a larger run, compare the complete p=13 CSV with the retained
twist-derived census. Independently compare the absolute-orbit weight
histogram and trace distribution with the complete PARI/GP p=13 census.
Keep any disagreements. At p=7, enumerate every field element for exact
point counts to audit the randomized counter's known assignment errors.

For each new completed prime, run 5,000 independent norm-one PARI/GP
point counts in a different field model. Check trace congruence, every
sample's positive-label membership, weighted totals, twist symmetry,
Hasse bounds, and conductor strata. These controls validate positive
membership and structural consistency; certify zero labels separately.

## Queue, resources, and failure handling

The authorized queue is
`53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199`.
Use four native census workers and sequential primes. This is a
correctness and class-label campaign on an ordinary host. Record wall
seconds for resource accounting, without a controlled timing ratio.

The sequential queue retains compressed results under
`/private/tmp/iso1-larger-primes-20261009`, where the available volume is
larger than the repository volume. Before each prime, require free space
of at least one GiB plus twice a conservative raw-CSV estimate of
`140*(p^3+1)` bytes. A storage-preflight refusal preserves the original
range as unfinished. After validation, compress with zstd, verify the
archive checksum, and compare its decompressed SHA-256 with the original
CSV before removing the queue's own uncompressed copy. Keep every receipt.

Record the exact command, source and executable hashes, field modulus,
seed law, thread count, raw stdout/stderr, exit code, CSV hash, and all
validation failures. Stop on a crash, disk exhaustion, or an explicit
recorded resource cap. Preserve partial output and its status. A launched
or queued prime is not a completed census. New ring-class hypotheses
will be registered with their training primes and exact formulas before
using larger completed primes as holdouts.

## Completed continuation evidence

The p = 53 census completed with 148,878 Hasse rows and 15,786,580
weighted representatives. Of 146,068 ordinary rows, 72,540 are weak;
0 of 73,034 depth-1 rows are weak, and 494 of 73,034 high-depth rows
are zero. All 5,000 independent GP controls (4,404 distinct absolute
traces) are in the positive set. The worker validated this completed
CSV and started p = 59; the remaining primes retain their queued status.
See `p53_absolute_validation.txt`, `larger_prime_source_freeze.txt`,
and the dated `larger_prime_status_20261008.tsv` snapshot. The live
status file remains under `/private/tmp/iso1-larger-primes-20261009/run2`.

## 2026-10-09 larger-field expansion and restart

The old temporary output directory and p59 worker were absent at the
start of this round. Preserve the old status as a historical snapshot.
The replacement finite launchd service uses a checked internal-volume
runtime and persistent receipts under
/Users/adamburan/Library/Application Support/crypto-iso1/census-20261009-run3.
It restarts at p59 and retains the complete requested prime list through
p199. See larger_fields_20261009/census_restart_status.tsv and
census_launchd.plist. The source supports an explicit relocated runtime
manifest; census arithmetic and the census executable are unchanged.

The user also requested larger primes and odd extension degrees. The
additional report verifies 17 paired point counts through log2(Q)=252
and total extension degrees 10 and 14, plus the prime-degree orbit
theorem. These controls do not complete the full trace census.
