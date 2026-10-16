# Order-1,201 code and evidence validation

- `rustfmt --check` passed for all touched Rust examples.
- `zsh -n run_order1201_v2.zsh` and `jq empty RUN_SPEC_V2.json`
  passed before child measurement.
- `cargo test --release --example` for the eight touched examples
  exited **0**; each example has zero embedded unit tests. The
  [complete command log](TOUCHED_EXAMPLES_TEST_V1.log) and
  [exit](TOUCHED_EXAMPLES_TEST_V1.exit) are retained. Functional
  correctness is supplied by the independent screen, full-fiber and
  84-child batch checkers, rather than by empty example test binaries.
- The selected full-fiber runner/checker both exited zero and agreed
  on the full map and 1,442,401 subgroup pair identities. The V2 batch
  and independent checker both exited zero; all 84 children, 1,785
  scalars and 12,486 relation witnesses passed.
- `actionlint .github/workflows/isogeny-fiber-four-order-replay.yml`
  passed after the historical source restoration step was added.
- `src/` has no changes relative to the preceding PR head
  `de340fa2161088a5188b841ec82278f5067f8aed`.

The required release library suite was run earlier in this PR; its
[validation receipt](../validation/VALIDATION_V1.md) retains the
4070-pass/13-failure/117-ignored full run and a later partial rerun,
plus passing focused socket and cover-test reruns. The failures were
sandbox socket denial and shared-host test timeouts outside this
change's paths. This remains a **red full-suite receipt**; touched
example compilation and independent arithmetic replay do not relabel
it green. A new source-library change or CI failure would require a
fresh affected validation.
