# Native execution commands

Working directory: the feature worktree rooted at
`/Volumes/SSD990/crypto/worktrees/codex-p192-p224-large-degree-isogenies-20261008`.
Source inputs and executable hashes are bound by `MANIFEST.json` after closeout.
The mathematical construction source is the new `isogeny_algos` library with
the exact-width dispatcher and P-224 characteristic-comparison repairs; the
supervisor does not change construction.

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

The first block records the initial search and partial replays. The corrected
phase rebuilt the CLI at default release optimization. The saved
`isogeny-algos-pre-p224-fix` identifies the earlier v2 trials; v3's monitor-race
receipt remains preserved. Before v4, the P-192 directory and both screening
files were copied byte for byte from v2. The supervisor checked each sealed
command and digest before reuse: 24 reused P-192 trials and 16 fresh P-224 trials.

```sh
cargo build --manifest-path isogeny_algos/Cargo.toml --locked --release --bin isogeny-algos
rustc --edition=2021 -C opt-level=3 isogeny_algos/examples/p192_p224_large_degree_search.rs \
  --extern isogeny_algos=isogeny_algos/target/release/deps/libisogeny_algos-a2c8eaec8fec19ed.rlib \
  -L dependency=isogeny_algos/target/release/deps \
  -o isogeny_algos/target/release/large-degree-search-monitor-fixed
isogeny_algos/target/release/large-degree-search-monitor-fixed \
  isogeny_algos/target/release/isogeny-algos \
  research/p192_p224_large_degree_isogenies_20261008/evidence-v4 \
  c70c32d486a3ac7531fe27f7193d9f09caa58344 --resume
CARGO_TARGET_DIR=research/p192_p224_large_degree_isogenies_20261008/target/algorithm-validation \
  cargo test --manifest-path isogeny_algos/Cargo.toml --locked --release --offline
cargo build --manifest-path research/p192_p224_large_degree_isogenies_20261008/Cargo.toml \
  --locked --release --offline --bins
research/p192_p224_large_degree_isogenies_20261008/target/release/replay \
  research/p192_p224_large_degree_isogenies_20261008/evidence-v4 docs/curves/registry.json
research/p192_p224_large_degree_isogenies_20261008/target/release/report \
  research/p192_p224_large_degree_isogenies_20261008 docs/curves/registry.json
research/p192_p224_large_degree_isogenies_20261008/target/release/catalogue-covers \
  docs/curves/registry.json docs/curves/covers.json docs/curves/cover-links.yaml
research/p192_p224_large_degree_isogenies_20261008/target/release/catalogue-covers \
  docs/curves/registry.json docs/curves/covers.json docs/curves/cover-links.yaml --check
research/p192_p224_large_degree_isogenies_20261008/target/release/catalogue-views
research/p192_p224_large_degree_isogenies_20261008/target/release/report \
  research/p192_p224_large_degree_isogenies_20261008 docs/curves/registry.json
pdftoppm -png -r 130 \
  research/p192_p224_large_degree_isogenies_20261008/SEARCH.pdf /private/tmp/p192-p224-final-page
```

`validation/algorithms-tests-final.log` has all 111 passing tests, including the
quadratic-wrapper check. The separate target kept the test build from replacing
the actual search executable. The example's two tests passed separately in
`validation/example-tests-final.log`. `validation/replay-final.log` verifies all
22 maps and raw hashes for all 40 distinct degree outcomes. P-192's additional
interrupted degree-149 invocation and earlier P-224 failures remain outside
those final degree counts. Executable digests and source revisions are assigned
to each evidence phase in `MANIFEST.json`. Closeout log paths are listed in
`EXECUTION_NOTES.md`; incomplete builds and failed launches remain preserved.

The dependent traits exporter was also attempted, and its existing root-library
compilation failure is preserved in `validation/curve-traits-build.log`:

```sh
cargo run --locked --release --offline --bin curve_traits -- build
research/p192_p224_large_degree_isogenies_20261008/target/release/report \
  research/p192_p224_large_degree_isogenies_20261008 docs/curves/registry.json --check-manifest
```

The last command is read-only and checks every byte length and SHA-256 digest
in the final manifest. The final manifest-generation invocation writes its
status to the terminal so it cannot change a log after hashing that log.
