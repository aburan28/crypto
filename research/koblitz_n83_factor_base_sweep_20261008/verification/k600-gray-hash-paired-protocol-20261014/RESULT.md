# Source-ready protocol for a matched K600 Gray/hash gate

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
vectors for each arm but **does not invoke them**. The
[paired supervisor](supervisor.py) requires a fresh exact 4,800-second grant
bound to the saved plan SHA256 and a fresh output path. It builds a separate
offline binary for each arm, launches one import and one search per arm under
the frozen Linux resource envelope, checks source/input/binary hashes between
stages, retains stage logs and container state, and conservatively charges an
entire 2,400-second arm once started. Producer failures stop the pair. The
[receipt audit](receipt_audit.py) independently replays the exact signed-orbit
difference and any candidate group equation, and checks receipt counters and
identities. Neither file invokes the worker on import. The plan's status is
`SOURCE_READY_AUTHORIZATION_PENDING`, and its selected-best field is null.

The grant schema is `n83.additional-local-budget/v7`, with purpose
`k600_gray_hash_paired_first_trial`, prior charge `23100`, additional amount
`4800`, the byte SHA256 of the saved plan, a fresh absolute `output_dir`, and
an explicit `user_approval_reference`. The archive directory must contain
the three pinned compressed objects named by their BLAKE3 digests. The run
requires a clean committed checkout and must be invoked explicitly using
`uv run --with blake3 python3 supervisor.py --authorization GRANT --output-dir
OUTPUT --archive-dir SCRATCH_DIR`. No such grant or run is recorded here.

With the three S3 archives downloaded under their compressed-BLAKE3 filenames
to `SCRATCH_DIR`, reproduce the source-only checks with:

```sh
uv run --with blake3 python3 research/koblitz_n83_factor_base_sweep_20261008/verification/k600-gray-hash-paired-protocol-20261014/preflight.py --archive-dir SCRATCH_DIR
uv run --with blake3 python3 -m unittest discover -s research/koblitz_n83_factor_base_sweep_20261008/verification/k600-gray-hash-paired-protocol-20261014 -p 'test_*.py' -v
```

The focused controls reject a modified parent audit, a stale or insufficient
grant, and malformed import receipts while pinning the matched base, distinct
residual sizes and search argument contract. This protocol adds no S3 objects
and does not alter the published factor-base archives.
