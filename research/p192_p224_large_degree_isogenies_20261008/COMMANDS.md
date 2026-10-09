# Native execution commands

Working directory: the feature worktree rooted at
`/Volumes/SSD990/crypto/worktrees/codex-p192-p224-large-degree-isogenies-20261008`.
Source inputs and executable hashes are bound by `MANIFEST.json` after closeout.
The mathematical construction source is the new `isogeny_algos` library with
the exact-width dispatcher repair; the supervisor does not change construction.

```sh
cargo build --manifest-path isogeny_algos/Cargo.toml --locked --release --bin isogeny-algos
rustc --edition=2021 -C opt-level=3 isogeny_algos/examples/p192_p224_large_degree_search.rs \
  --extern isogeny_algos=isogeny_algos/target/release/deps/libisogeny_algos-a2c8eaec8fec19ed.rlib \
  -L dependency=isogeny_algos/target/release/deps \
  -o isogeny_algos/target/release/large-degree-search-resumable
isogeny_algos/target/release/large-degree-search-resumable \
  isogeny_algos/target/release/isogeny-algos \
  research/p192_p224_large_degree_isogenies_20261008/evidence-v2 \
  c70c32d486a3ac7531fe27f7193d9f09caa58344 --resume
cargo build --manifest-path research/p192_p224_large_degree_isogenies_20261008/Cargo.toml \
  --locked --release --offline --bins
research/p192_p224_large_degree_isogenies_20261008/target/release/replay \
  research/p192_p224_large_degree_isogenies_20261008/evidence-v2 \
  docs/curves/registry.json p192:73
research/p192_p224_large_degree_isogenies_20261008/target/release/replay \
  research/p192_p224_large_degree_isogenies_20261008/evidence-v2 \
  docs/curves/registry.json p192:193
cargo test --manifest-path isogeny_algos/Cargo.toml --locked --release
cargo test --release --lib
```

The direct supervisor build avoided Cargo's concurrent test-artifact lock.
For portable reproduction use Cargo's named example instead; the exact rlib
filename above is the observed local artifact. The first successful launch used
`large-degree-search-supervised` without `--resume`, before resume support was
added. Completed child commands and public curve parameters are recorded in
each receipt; resumed receipts are reused only after command and byte-digest
checks. Neither launch uses the abandoned alternate optimization-level-1 build.

After the complete search finishes, run the replay without its final selection
argument, generate the report/target registration, regenerate covers and joins,
and then regenerate the manifest. Those closeout commands remain pending until
their logs and final artifacts are present.
