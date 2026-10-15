# Source-only protocol for a matched K600 Gray/hash gate

The [frozen plan](plan.json) names two exact primary-curve N83 arms. Both use
the same 600-column `public_x_sequential` small base (seed `2026100802`),
the same public fixture-zero point, first-trial coefficients, 3+2 equation,
5,000-conflict cap, 32-model cap, 1,100-second worker watchdog, and
one-CPU/four-GiB/no-swap search envelope. The only factor-base input changed
between arms is the K600 pool: `public_x_gray_prefix` versus `public_x_hash`,
both seed `2026100801`. The archives, S3 URIs and content digests are pinned
by the [live S3 audit](../pilot-pair-live-s3-20261014/RESULT.md).

| Pool | Residual signed orbits | Source total variables | Source domain clauses | Conditional 50%-allowance floor for 600 rows |
| --- | ---: | ---: | ---: | ---: |
| Gray-prefix | 97 | 95,068 | 183,314 | 33,975 queries |
| Hash | 600 | 97,080 | 269,830 | 889 queries |

The capacity figures come from the pinned existing worker formula for the
five-summand 3+2 shape. The query floors are necessary counting bounds for
uniform nonidentity queries and a specified 600-row goal, not measured hit
rates, SAT throughput or sufficient work estimates. One bounded search call
per arm could test solver access; even two successful calls would not measure
relation yield, useful rank, linear algebra, a fresh individual logarithm, or
complete cold runtime. The pool policies also differ, so an observed timing
gap would describe these exact configurations rather than prove that orbit
overlap alone caused it.

The plan models two independently charged cold arms at 2,400 seconds each:
600 seconds build, 360 import, 1,200 search and 240 cleanup per arm. This is
**4,800 additional seconds requested, not authorized**. The prior local N83
search ledger is 23,100 seconds; 27,900 is a hypothetical cumulative charge
only if a later grant is used in full. The earlier 2,400-second grant has
already been fully charged. No N83 solver search, relation, rank or runtime
measurement was produced by this source-only preflight.

The [preflight](preflight.py) checks byte-pinned manifest, corpus, S3 audit,
support screen, capacity source, worker source and Cargo lock, replays all
three local archive bytes through the live audit, recomputes both capacity
rows and rank-query floors, and compares the generated JSON to the saved plan.
It emits exact worker `verify-cross-pilot` and `run-cross-pilot` argument
vectors for each arm but **does not invoke them**. A new pair-specific
supervisor, independent result audit, clean-source/binary receipts, explicit
authorization and stopping rules are still required before the plan is
executable. The plan's status remains
`SOURCE_READY_SUPERVISOR_AND_AUTHORIZATION_PENDING`, and its selected-best
field is null.

With the three S3 archives downloaded under their compressed-BLAKE3 filenames
to `SCRATCH_DIR`, reproduce the source-only checks with:

```sh
uv run --with blake3 python3 research/koblitz_n83_factor_base_sweep_20261008/verification/k600-gray-hash-paired-protocol-20261014/preflight.py --archive-dir SCRATCH_DIR
python3 -m unittest discover -s research/koblitz_n83_factor_base_sweep_20261008/verification/k600-gray-hash-paired-protocol-20261014 -p 'test_*.py' -v
```

The focused controls reject a modified parent audit and pin the matched
base, distinct residual sizes and search argument contract. This protocol
adds no S3 objects and does not alter the published factor-base archives.
