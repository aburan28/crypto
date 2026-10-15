# First retained-N83 mixed-orbit producer gate

This protocol is **prepared, not authorized or launched**. The prior local
allowance is fully charged at 11,700 seconds. The supervisor requires a new,
run-specific `n83.additional-local-budget/v1` record naming the exact fresh
build and run directories, purpose `mixed_orbit_two_residual_first_gate`,
`prior_pilot_charged_seconds` `11700`, at least **2,400 additional seconds**,
and a user approval reference. A source file or synthetic test record is not
an authorization. Do not run either command without that separate grant.

The frozen first arm is primary Koblitz `a=0`, seed `2026100801`,
`public_x_hash`: K=1,182 as the small signed-Frobenius base and the exact
K=2,048 envelope, leaving 866 residual signed orbits. Use the existing
fixture-zero public target whose canonical corpus BLAKE3 is
`831eaa7be2387634fd610f6706103b2af9b4e274f0db59923c0280af3cbed74b`.
Search one deterministic public equation with three small and two residual
summands, ordered first small pair, at most 150,000 variables, 12,000,000
domain clauses, 32 models and 5,000 conflicts. These settings are a capacity
and partial-yield gate; they do not establish a factor-base runtime ranking.
The independently derived first-trial coefficients are
`a=1833078656608906180176350` and `b=163089659398503559150141` modulo
the pinned subgroup order. The outer event audit rejects a different source
equation before considering any partial or graph row.

The frozen local objects must be downloaded from their retained S3 URIs and
checked against these SHA-256 digests before launch:

| Role | S3 object | SHA-256 |
| --- | --- | --- |
| Small | `s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/v2-size-frontier/a0/objects/76e1cbe4ea5a9d26e5b24a2d27f9401271dc471846529190a0416d8e09837a92.jsonl.gz` | `eb342b574e244e1e67e8e78c805f33ee88edcdb4f3a14d4114e6f413beede380` |
| Envelope | `s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/v2-size-frontier/a0/objects/82990cda5d57b1c4c4b635673a881ad233c24f707e7ac644a11478aeddddf851.jsonl.gz` | `aacb62784dee0c9f5e06328dbb85a4b6a1c7243e778328da82900e7284f3a74a` |

Build only from a clean commit with the pinned offline Linux arm64 image
`sha256:7c4ae649a84014c467d79319bbf17ce2632ae8b8be123ac2fb2ea5be46823f31`.
The build gets four CPUs, 8 GiB, no network and a 450-second outer wall; its
conservative charge is 600 seconds. The `verify` stage gets one CPU, 4 GiB,
zero swap, no network and a 360-second outer wall. Only after a valid import
receipt and confirmed container cleanup may the search start. Search gets the
same resource isolation, a 1,200-second outer wall and a 1,150-second worker
watchdog. A 240-second reserve covers cleanup, frozen-input checks, event
replay and receipts. Full planned conservative charge is 2,400 seconds;
only launched stages and the applicable cleanup reserve are charged in the
budget receipt. A timeout, OOM, invalid model, failed import or malformed
receipt remains visible and does not trigger an automatic retry.

The entry point is
`research/koblitz_n83_factor_base_sweep_20261008/mixed_orbit_supervisor.py`.
Its `build` mode requires `--authorization`, `--build-dir` and `--output-dir`;
`run` additionally requires `--small-object` and `--envelope-object`. Both
directories must be absolute paths outside the checkout, with the build
directory fresh for `build` and the output directory fresh for `run`. The
same authorization file binds both stages. Preserve its bytes, the committed
source, Linux ELF, image ID, exact commands, raw stdout/stderr, container
state and cleanup, worker config/result/events, independent group replay,
outer receipt and budget receipt. Charge download and artifact handling in
the eventual complete cold accounting even when this gate only reports stage
metrics.

Before using any relation-yield result to choose a base, repeat comparable
gates for competing bases and policies, then measure graph completion,
independent rank growth, sparse linear algebra, one individual logarithm,
complete cold single-target runtime and matched rho under calibrated,
isolated resources. Censored gates cannot rank total runtime.
