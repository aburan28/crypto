# Explicit residual-pool screen for the retained N83 bases

The frozen v1 grid names zero, one or two large primes, but that count alone
does not define the residual-point domain. This screen fixes one finite policy
using objects already stored in S3: for the same curve arm, public-x policy
and seed, take the first K=64 or 256 signed-Frobenius orbit columns as the
small base, and take the later columns through K=256 or 600 as the residual
pool. This is an exploratory pool policy for a later design version. It does
not change the frozen v1 grid or turn its large-prime labels into executed
solver arms.

The bounded checker read the first JSON line of all 54 locally retained
compressed objects. For all 18 curve/policy/seed groups it verified that
K=64 representatives and accepted candidate indices prefix K=256, which
prefix K=600. Each object header was bound to its manifest row, and the
existing 54-object generic replay and S3 round-trip receipts were required
to pass. The resulting residual pools have exactly 166 times the difference
in orbit columns as distinct point records:

| Small K | Envelope K | Small points B | Residual points L |
| ---: | ---: | ---: | ---: |
| 64 | 256 | 10,624 | 31,872 |
| 64 | 600 | 10,624 | 88,976 |
| 256 | 600 | 42,496 | 57,104 |

For m total summands, exactly j residual points and repetitions allowed,
the number of unordered point multisets is

    C(B + m - j - 1, m - j) * C(L + j - 1, j).

The cumulative at-most-two count adds j=0,1,2. Every multiset has one fixed
group sum. For a uniform nonidentity target in a prime-order subgroup of
order r, the expected number of candidate decompositions is at most the
cumulative count divided by r-1; Markov gives the same expression clipped
at one as a rigorous hit-probability ceiling. Identity-summing multisets
only make that ceiling looser. This does not measure the solver's ability
to find a partial or the graph's ability to close a cycle.

Primary a=0 ceilings for at most two residuals are:

| Small K → envelope K | m=4 | m=5 | m=6 |
| --- | ---: | ---: | ---: |
| 64 → 256 | 1.471146e-8 | 4.946403e-5 | 0.1272816 |
| 64 → 600 | 9.997757e-8 | 3.472970e-4 | 0.9118867 |
| 256 → 600 | 9.672326e-7 | 0.01231350 | 1 (vacuous) |

The machine receipt retains exact integers and all 90 cases for both curve
arms, three pool pairs, five arities and three residual counts. Nine stored
variants per arm share each pair's B and L but are dependent construction
variants, not nine independent yield observations. The diagnostic a=1 arm
remains secondary.

An independent receipt checker recomputes every integer case and verifies
the manifest, generic-replay and S3 round-trip bindings and the wall-budget
digest. Its passing log is verification/lp-pool-independent-receipt-check.log.
It reads receipts, not compressed point objects; the bounded producer is the
only run in this update that re-opened all 54 headers.

The checker used a six-second child wall cap. It finished in 0.419 seconds;
the budget receipt charges a further one-second process-overhead allowance,
bringing the conservative active pilot total to 3,593.093 of 3,600 seconds.
The receipt is a support screen, not a timed relation, matrix or DLP run.

A cold comparison must charge construction or retrieval of both the small
base and its envelope, residual-pool indexing, every unsuccessful partial,
graph filtering and rank, linear algebra, individual logarithm and
verification. A source-pinned producer for this exact pool and an independent
replay are still needed. No factor-base runtime winner follows from these
ceilings.

The bounded command refuses to overwrite its immutable outputs. A new
retained-data replay needs a separately versioned study path and receipt.
The pure combinatorial and prefix controls run with:

    python3 -m unittest discover -s research/koblitz_n83_factor_base_sweep_20261008 -p test_large_prime_pool_screen.py -v

The retained artifact is pilot-01/large-prime-pool-screen.json, paired with
pilot-01/large-prime-pool-budget.json. The exact source is
large_prime_pool_screen.py; the independent checker is
verify_large_prime_pool_screen.py.
