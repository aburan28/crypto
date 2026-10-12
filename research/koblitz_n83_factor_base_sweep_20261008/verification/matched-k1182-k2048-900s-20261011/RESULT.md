# Matched K=1,182 and K=2,048 bounded N83 search

Both retained primary signed-Frobenius bases reached `UNKNOWN_solver_cap` at
exactly 5,000 conflicts on the **same public point**. Each one-trial,
five-summand search produced zero SAT models, zero natural relations and rank
zero. This is a pair of censored stage observations. Neither factor base has
a measured complete cold index-calculus runtime, and no minimum-runtime base
is selected.

The source commit was `831cedf9eea92ab21751965b9769687855bb898a`; the
single Linux worker had SHA-256
`a195936edd6ed15fa96201dd69b54f62c67e99871e3ba7d02c2f12713c753c7d`.
Both arms used public fixture zero, curve `a=0`, policy `public_x_hash`, seed
`2026100801`, the unrestricted Frobenius-orbit-phase domain, five summands,
one trial, 150,000 variable and 1,000,000 domain-clause caps, at most 32
models, 5,000 conflicts, one CPU, 4 GiB hard memory, zero swap, no network,
and a 280-second outer wall cap. The worker's target point in both receipts
was `(6261116474805304534310590, 4505575220639647215230954)`. The
supervisor's original-curve group replay passed for each query.

| Base | Stored points | Outer wall, L0 | Base import, L0 | Worker solver stage, L0 | Cgroup peak | Conflicts | Models | Relations | Rank |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| K=1,182 | 196,212 | 188.821 s | 2.264 s | 185.709 s | 397,496,320 B | 5,000 | 0 | 0 | 0 |
| K=2,048 | 339,968 | 168.456 s | 4.310 s | 163.542 s | 556,707,840 B | 5,000 | 0 | 0 | 0 |

The K=2,048 outer wall was 20.364 seconds shorter in this one shared-host
observation. Neither arm emitted a relation, so that difference is not a
relation-yield or total-runtime advantage. Column logs, sparse linear
algebra at full subgroup width, individual logarithm and matched rho were not
run in this gate. The exact stored-object comparison reconfirmed that all
1,182 representatives and 196,212 point rows of the smaller base are a
prefix of the larger; its suffix has 866 orbits and 143,756 point rows. This
is object geometry, not a measured large-prime yield.

The K=1,182 and K=2,048 compressed objects are at
`s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/v2-size-frontier/a0/objects/76e1cbe4ea5a9d26e5b24a2d27f9401271dc471846529190a0416d8e09837a92.jsonl.gz`
and
`s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/v2-size-frontier/a0/objects/82990cda5d57b1c4c4b635673a881ad233c24f707e7ac644a11478aeddddf851.jsonl.gz`.
Their previously retained generic point replay and upload/download receipts
both passed. This gate used local copies with independently checked compressed
SHA-256 values `eb342b574e244e1e67e8e78c805f33ee88edcdb4f3a14d4114e6f413beede380`
and `aacb62784dee0c9f5e06328dbb85a4b6a1c7243e778328da82900e7284f3a74a`.
It did not upload a new object. The replay uses a different arithmetic
backend in the same repository and host; external independent replay remains
open.

The Linux build passed in 136.034 seconds of process wall with its source and
cleanup checks. The two search supervisors accepted the worker, source,
object and fixture hashes, replayed the query, and verified their containers
were removed. The independent [audit](audit.py) reruns receipt and hash
checks, complete object-prefix traversal, original-curve query replay and
container-absence checks. Its [machine result](evidence-01/audit.json) is
`PASS_bounded_receipts_only`; [raw evidence](evidence-01/) includes build,
worker, outer, panel, exact outer commands, the frozen source archive and
evidence-manifest receipts. The five source files attested by the worker are
retained in `evidence-01/source-snapshot.tar.gz`; the audit checks their exact
hashes and, when the original source commit is locally available, its Git
blobs. That source commit remains reachable on branch
`codex/n83-ordered-first-pair-20261011`. The Linux worker and compressed
objects remain outside Git at their frozen local/S3 locations.

Validation: the original worker source's study Python suite passed 39/39, and
its required full `cargo test --release --lib` recheck passed 3,898 tests with
110 ignored. An initial sandboxed full run had 3,889 passes and nine unrelated
local-socket tests blocked by `Operation not permitted`; the loopback-enabled
full recheck passed them. After this evidence was ported onto current main
`c62555b7a1746490dba707589d791669a20bb0d0`, the study suite again passed
39/39 and the full release library suite passed 3,986 tests with 112 ignored.
All logs and the [validation receipt](evidence-01/validation/validation.json)
are retained. The receipt audit passed with read-only Docker access both with
the original Git object and with an empty Git database using only the source
archive. The source/document diff check was clean; raw logs retain their exact
output.

The predeclared [protocol](PROTOCOL.md) charges 280 seconds for each launched
build/search stage and 60 seconds for cleanup and receipts, regardless of
their shorter actual walls. The full **900-second** extension is therefore
charged, increasing the cumulative conservative local charge from 10,800 to
**11,700 of 11,700 seconds**. No further retained N83 experiment is authorized
by this grant. The next discriminating gate needs measured natural relation
yield and rank novelty under a useful solver budget, followed by complete
precompute, column-log and individual-log accounting against matched rho.
