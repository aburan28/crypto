# N83 signed-orbit size frontier: exact design extension

The frozen v1 design stops at 900 public-x orbit columns. The primary
`a=0` subgroup has order `2417851639230796216685689`. For a successfully
constructed signed-Frobenius base with `K` full, distinct orbit columns,
`B=166K` point records. The exact number of unordered `m`-point multisets
with repetition is `M=binomial(B+m-1,m)`. For a uniformly chosen nonidentity
target, `min(1,M/(r-1))` is an **upper bound** on the probability of any
full-smooth decomposition. A value of one makes the bound uninformative; it
does not establish a relation, solver throughput, rank gain or runtime.

The existing support receipt proves that this bound first reaches one at
`K=1,182` for `m=5` and `K=16,627` for `m=4` on the primary arm. The new
sizes `1,182, 2,048, 4,096, 8,192, 16,627` bracket those transitions.
The table includes the prior 600 and 900 sizes for scale. Every integer and
fraction is in `size-frontier-v2.json`; decimal ceilings are rounded.

| K | Conditional points B | m=4 hit ceiling | m=5 hit ceiling | Full pair multisets |
| ---: | ---: | ---: | ---: | ---: |
| 600 | 99,600 | 1.695987e-6 | 3.378543e-2 | 4,960,129,800 |
| 900 | 149,400 | 8.585764e-6 | 2.565495e-1 | 11,160,254,700 |
| 1,182 | 196,212 | 2.554316e-5 | 1 | 19,249,672,578 |
| 2,048 | 339,968 | 2.302072e-4 | 1 | 57,789,290,496 |
| 4,096 | 679,936 | 3.683283e-3 | 1 | 231,156,822,016 |
| 8,192 | 1,359,872 | 5.893227e-2 | 1 | 924,626,608,128 |
| 16,627 | 2,760,082 | 1 | 1 | 3,809,027,703,403 |

The pair column counts multisets a materialized complete pair index would
have to address; it is an operation/index-size diagnostic, not a lower bound
for every possible splitting, SAT, FES or sparse algorithm. For the more
restrictive model of uniform queries yielding at most one full-smooth row,
the exact Markov rank ceiling cannot reach one half before 23,137,309
four-summand queries at `K=1,182`, 69,504 at `K=8,192`, and 16,627 at
`K=16,627`. Those are necessary query counts only. Fixed public targets,
dependencies among rows, solver success and large-prime cycles are outside
this bound.

The versioned design appends 270 new public-x base specifications across two
curve arms, three policies, five sizes, three seeds and three closure choices.
It retains every v1 base and solver ordinal as a prefix, then pairs each new
base with all 43,200 frozen solver-axis combinations. This makes 57,024,000
addressable design tuples, of which 11,664,000 are new and unexecuted. An
ordered source-availability rule assigns exactly one disposition to every
tuple. It identifies the branch's wide large-prime row adapter while marking
its missing partial producer; a disposition never means an arm ran. The
`--case` command decodes any ordinal without generating all tuples:

    python3 research/koblitz_n83_factor_base_sweep_20261008/size_frontier_v2.py --case 45360000

The receipt binds the v1 design, the support receipt and its panel manifest by
SHA-256. The generator refuses to overwrite changed output. Its five tests
check exact threshold minimality, v1 ordinal-prefix preservation, bijective
mixed-radix solver decoding, new-size addressability, exact support fractions,
the one-row rank-query arithmetic, and aggregate dispositions against the
frozen Rust counts:

    python3 -m unittest discover -s research/koblitz_n83_factor_base_sweep_20261008 -p test_size_frontier_v2.py -v

On the current branch, the release library suite passed 2,249 tests (94
ignored), the N83 exporter example passed 16, the study Python suite
passed 13 and the boundary Python suite passed 16. The v2 generator was rerun
against its saved JSON without a byte change. These are design and source
checks, not N83 performance measurements.

This extension constructs **no new factor-base objects** and makes no
single-target cold timing claim. Successful full-orbit construction at larger
K, a capacity check of the current buffered exporter or a streaming output
path, source-pinned relation production, natural rank and fully charged cold
comparison remain required. The declared
S3 destination for future versioned objects is
`s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/v2-size-frontier/`.

The exporter now has a one-object path for the five added sizes. It accepts
only the two pinned curve arms, the three public-x policies, the three frozen
seeds and signed-Frobenius closure. Construction requires a clean committed
worktree and checks that the compiled exporter matches its on-disk source.
A fresh directory receives a
content-addressed compressed object and a manifest bound to this v2 design;
a process-wall watchdog records `UNKNOWN_budget` if construction exceeds its
cap. The separate bounded replay uses generic multi-limb curve arithmetic to
check source points, cofactor projection, subgroup and Frobenius identities,
every point label, closure and the point-set hash. `upload` requires that
replay, uploads to the versioned S3 prefix, downloads the object to rehash its
bytes, and uploads the manifest and receipt. The commands, after a new
construction/replay budget and a memory guard are authorized, are:

    cargo run --release --example koblitz_n83_factor_base_export -- v2-construct-one NEW_DIRECTORY 0 public_x_hash 1182 2026100801 BUDGET_SECONDS
    cargo run --release --example koblitz_n83_factor_base_export -- v2-replay-one NEW_DIRECTORY BUDGET_SECONDS
    cargo run --release --example koblitz_n83_factor_base_export -- upload NEW_DIRECTORY

The size and destination gate and the v2 replay schema pass the existing
small public fixture; out-of-grid publication fails before any AWS call. The
CLI also rejects a valid v2 construction request from an unfrozen worktree
before creating its output directory, and rejects a v1 manifest at the bounded
v2 replay entry point. These negative controls do not exercise a larger base.
Larger-object construction, its memory use, generic replay time and S3
round-trip have **not** been run under the exhausted pilot budget. The
watchdog bounds process wall only, so a future run also needs an external
memory limit and a fresh run directory.
