# N83 mixed-orbit worker source gate

The worker now connects the two pinned K=1,182/K=2,048 size-frontier archives
to exact residual-suffix construction, bounded one/two-residual chained-S3
search, and full-width large-prime graph replay. The `verify` command checks
both complete archives and their prefix relationship. The `run` command
requires the matched public-corpus hash, a 4 GiB zero-swap Linux cgroup, explicit wall/trial/conflict/model
caps, and a fresh output directory. It records every public source equation,
search classification, accepted partial and graph outcome. An invalid SAT
model is a producer failure; preflight, search and graph caps remain unknown.

This is a **source-only** gate. Neither `verify` nor `run` was invoked on the
retained N83 archives. There are no new partial-relation, graph-completion,
rank, linear-algebra or complete cold-runtime measurements. The source tests
check pinned-archive rejection and deterministic public trial coefficients;
they do not validate a successful N83 import or search.

The example's release build and three focused tests passed. The full release
library suite passed 4,000 tests with zero failures and 113 ignored. The study
Python suite passed 39/39. The touched Rust source passed `rustfmt --check`,
and the diff passed `git diff --check`. The [machine receipt](validation.json)
binds the source SHA-256 and deterministic compressed raw logs.
