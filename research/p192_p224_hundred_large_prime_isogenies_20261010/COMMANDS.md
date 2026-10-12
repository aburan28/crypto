# Native campaign commands and accounting

Working directory: the isolated PR #1606 worktree. Source protocol commit
10743ca59 precedes constructor execution. Producer/verifier source commit
49a280dc60972c9d16b2205b4652cc6c5d87284e was clean when the binary was built.

The existing dedicated CI binary performed both complete prime screens:

```sh
isogeny screen --curve p192 --from 1010 --to 1000000 --max-order 8
isogeny screen --curve p224 --from 1010 --to 1000000 --max-order 8
isogeny screen --curve p192 --from 1010 --to 1000000 --max-order 64
isogeny screen --curve p224 --from 1010 --to 1000000 --max-order 64
```

Maintainer validation and build (runtime users execute the resulting binary):

```sh
cargo test --release --manifest-path tools/isogeny-cli/Cargo.toml
cargo test --release --manifest-path isogeny_algos/Cargo.toml
cargo clippy --manifest-path isogeny_algos/Cargo.toml --locked --all-targets -- -D warnings
cargo test --release --lib
cargo clippy --release --manifest-path tools/isogeny-cli/Cargo.toml --bin isogeny -- -D warnings
cargo build --release --locked --manifest-path tools/isogeny-cli/Cargo.toml --bin isogeny
```

The root library and additional strict CLI lint outcomes remain retained failures.
Research outputs and build caches use separate target directories. Both actual
dedicated executable invocations ran concurrently, with fresh output directories:

```sh
isogeny walk --curve p192 --ell 2341 --count 50 --out evidence/p192
isogeny walk --curve p224 --ell 1117 --count 50 --out evidence/p224
```

The exact commands for every worker, input model hashes, stdout/stderr hashes,
exit status, sampled peak RSS, and constructor elapsed intervals are in each
step's search.json and ell-N.receipt.json. Walk manifests additionally record
independent replay and complete step intervals. Two walks share the host;
summed intervals are aggregate work intervals, not campaign wall time or a
throughput comparison. Independent replay sits outside the constructor's
timeout and RSS cap. Map counts exclude replay duplication.

Native walk-controls copies the first two sealed P224 steps and checks a valid
two-step replay, a modified next-source model with its manifest hash updated,
and a modified transported generator with its manifest hash updated. The latter
two must fail mathematical/source binding checks despite their updated hashes.
Its log, inputs and stdout/stderr outputs are retained as separate controls and
never counted as newly found maps.

The report builder runs offline from its native locked package, registers only
verified destinations, writes the diagram/PDF and seals the artifact manifest.
An initial report build failed DNS resolution; its log remains preserved before
the successful offline retry. Native catalogue-covers, catalogue-views and
canonical-roster refresh and validate the canonical generated copies. The
original studies and root IC measurement bytes are preserved.
