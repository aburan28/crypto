# Pinned K=600 mixed-orbit worker import

The existing N83 mixed-orbit worker now accepts `verify-alt` and `run-alt` for
the K=600/K=2,048 nested pair. The original `verify`, `run`, `verify-control`
and `control-run` retain the K=1,182/K=2,048 pair. Both natural search modes
use the same public corpus, trial coefficients, S3 solver path, first-pair
symmetry, bounded search settings, graph eliminator and result schema. The
selected small-object and manifest hashes are written to each run's config
and result. The control modes remain bound to their original known-answer
fixture.

The K=600 object is selected by its exact compressed BLAKE3 digest from the
original 54-row pilot panel. The importer checks the pinned manifest bytes,
panel schema and count, unique selected row, S3 URI, compressed and plain
digests, every orbit-labelled point, subgroup and Frobenius identities, and
the exact residual suffix reconstructed from the K=2,048 envelope. `run-alt`
retains the existing 4 GiB, zero-swap Linux cgroup gate and requires an
external budget supervisor. The alternative mode does not fetch an object.

The following imports succeeded using local copies of the already stored S3
objects. Neither invocation called the relation solver:

The [machine import receipt](import-receipts.json) records the exact worker
source SHA-256, manifest and compressed-object SHA-256 digests, and the raw
JSON response from each importer invocation.

| Mode | Small BLAKE3 | Small points/orbits | Residual points/orbits | Status |
| --- | --- | ---: | ---: | --- |
| `verify-alt` | `b4dda45e6c3edfca673db6f3f82c16cce8c860f5e62c092155258df399a0363e` | 99,600 / 600 | 240,368 / 1,448 | `PASS`, `solver_search_executed=false` |
| `verify` | `76e1cbe4ea5a9d26e5b24a2d27f9401271dc471846529190a0416d8e09837a92` | 196,212 / 1,182 | 143,756 / 866 | `PASS`, `solver_search_executed=false` |

The pilot K=600 manifest SHA-256 is
`5bd862987844f709e03beeab9fd265e3d2024011bcadcb82f5af9c72018df6a6`;
its local compressed object SHA-256 is
`3e85c9e419fcddd3c9faed3c84470bf96fcad9178e48e1d73629369335841ecf`.
The K=1,182 manifest and object SHA-256 digests are respectively
`2735f6960233749ac887f68705fd514dfc003e2212d722a8e50e5e9738aadcfa`
and `eb342b574e244e1e67e8e78c805f33ee88edcdb4f3a14d4114e6f413beede380`.
The shared K=2,048 manifest and object SHA-256 digests are respectively
`68b43524198394a2f0d112f86df0d1c1f69ea1ff817695cf21efe2767b7d17b8`
and `aacb62784dee0c9f5e06328dbb85a4b6a1c7243e778328da82900e7284f3a74a`.
The K=600 S3 HEAD and upload receipt are retained in the adjacent
[prefix screen](../k600-k2048-prefix-screen-20261012/RESULT.md).

For the pinned three-small/two-residual shape, the source formula in
`frobenius_orbit_domain_size(83, m, representatives)` gives the following
finite-domain preflight totals. These are arithmetic counts, not measured SAT
resource use:

| Pair | Domain variables | Domain clauses | XOR rows |
| --- | ---: | ---: | ---: |
| K=600/K=2,048 | 15,647 | 415,686 | 3,320 |
| K=1,182/K=2,048 | 16,811 | 465,738 | 3,320 |

The 465,738-clause count agrees with the existing K=1,182 search preflight.
The K=600 path still needs a measured solver attempt before its lower clause
count can be interpreted as a runtime effect.

Validation: the touched release example passed four tests; the study Python
suite passed 46 tests and the boundary Python suite passed 43; both exact
Rust import modes passed and their raw responses are in the machine receipt.
The full release library suite passed 4,000 tests with 113 ignored when
localhost binding was permitted. Its initial network-restricted run had nine
unrelated localhost-test failures at `TcpListener::bind` with
`PermissionDenied`; the identical full command passed with that sandbox
restriction removed. `git diff --check` passed.

This is source and import evidence only. It adds no K=600 relation-yield,
graph-cycle, rank, sparse-linear-algebra or complete-runtime measurement. In
particular, smaller columns alone cannot rank the two pairs by total index
calculus runtime. A future natural-target comparison needs a separately
authorized bounded run on an isolated host, followed by replay of the
config, event stream, result and complete phase accounting.
