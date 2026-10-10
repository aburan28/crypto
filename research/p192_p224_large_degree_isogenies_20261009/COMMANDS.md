# Native continuation execution

Working directory is the isolated codex-p192-p224-large-degree-isogenies-20261008
worktree. The original study remains frozen; this study is a distinct round.

The full standalone release suite, including three new exact modular-series
checks, passed 114 tests. The native supervisor's two structural/order checks
passed. Extension arithmetic passed four tests covering reducibility rejection,
field inverses and Frobenius, the source-sized fields, and an exhaustive
small-field point count independent of the standard trace recurrence.

```sh
CARGO_TARGET_DIR=/private/tmp/p192-p224-continuation.f9sck4/algorithm-validation \
  cargo test --manifest-path isogeny_algos/Cargo.toml --locked --release --offline
```

Full-Hecke construction executable:
`/private/tmp/p192-p224-continuation.f9sck4/bin/isogeny-algos-sieve-97c80540b`.
Its source algorithm revision is 97c80540b. Four serial supervised invocations
use seed 1 and the published source orders: p192/509/600s, p224/521/600s,
p192/1021/1800s and p224/1031/1800s. `ISOGENY_STAGE_TRACE=1` records stage
completion on stderr. Each run has its own directory and `search.json`.

The kernel-first helper source was frozen at 2fe16c12b. It was built directly
with rustc at optimization level 3 against the tested isogeny_algos rlib. The
helper is `/private/tmp/p192-p224-continuation.f9sck4/kernel-first-v3`.
`supervisor-single` invokes it with the same CLI options and --single-map;
the latter requires exactly one map and binds KERNEL_PROTOCOL.md. The original
full-Hecke supervisor continues to require exactly two maps. Kernel-first
cases are p224/1471/1800s and p192/10453/1800s, run serially after all four
full-Hecke outcomes have completed.

The independent per-map verifier is built as follows:

```sh
CARGO_TARGET_DIR=/private/tmp/p192-p224-continuation.f9sck4/replay \
  cargo build --manifest-path research/p192_p224_large_degree_isogenies_20261009/Cargo.toml \
  --locked --release --offline --bin replay-one
/private/tmp/p192-p224-continuation.f9sck4/replay/release/replay-one \
  KERNEL_EVIDENCE_DIR docs/curves/registry.json
```

The verifier retains all mathematical kernel, codomain, exact map-identity and
public subgroup transport checks from the first round. Only the explicitly
declared enumeration scope differs: it requires one map, one specified
eigenline and degree coverage PARTIAL. It cannot report complete two-eigenline
enumeration. All failed launches/builds and bounded outcomes are preserved.

High-degree exact verification additionally uses validated reciprocal
polynomial reduction, memoized division recurrences and block modular
composition. Twelve tests compare these with the unchanged reference routines.
The initial eager degree-10453 replay exited 143 before a certificate; the
separately preserved memoized replay completed all checks. The exact fixture
control was replayed again before target use.

```sh
CARGO_TARGET_DIR=/private/tmp/p192-p224-continuation.f9sck4/replay \
  cargo test --manifest-path research/p192_p224_large_degree_isogenies_20261009/Cargo.toml \
  --locked --release --offline --bin verification-poly-checks --bin verification-kernel-checks
/private/tmp/p192-p224-continuation.f9sck4/replay/release/continuation-report \
  research/p192_p224_large_degree_isogenies_20261009 docs/curves/registry.json
research/p192_p224_large_degree_isogenies_20261008/target/release/catalogue-covers \
  docs/curves/registry.json docs/curves/covers.json docs/curves/cover-links.yaml
research/p192_p224_large_degree_isogenies_20261008/target/release/catalogue-covers \
  docs/curves/registry.json docs/curves/covers.json docs/curves/cover-links.yaml --check
research/p192_p224_large_degree_isogenies_20261008/target/release/catalogue-views
/private/tmp/p192-p224-continuation.f9sck4/replay/release/continuation-report \
  research/p192_p224_large_degree_isogenies_20261009 docs/curves/registry.json
/private/tmp/p192-p224-continuation.f9sck4/replay/release/continuation-report \
  research/p192_p224_large_degree_isogenies_20261009 docs/curves/registry.json --check-manifest
```

The final report/manifest invocation prints to the terminal so no study log
changes after being hashed. Canonical catalogue outputs and the old traits
file are bound, along with source files, raw evidence and executed binaries.
The root release-library check remains a separately recorded failed gate.
