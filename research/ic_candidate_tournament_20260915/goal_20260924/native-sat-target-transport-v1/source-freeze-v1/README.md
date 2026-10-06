# Target-free SAT source/build freeze

This input fixes a bounded n17 external-CryptoMiniSat target policy before a
fresh public point exists. The exact 512-query ordinary preparation and the
nine-role marked-CMS control are prerequisites. The preparation paths identify
original local evidence, so another machine must verify those same sealed
bytes at its own canonical paths and issue a new registration.

`icprog sat-target-freeze` copies the committed Rust source, vendors locked
offline dependencies, builds the worker, controller and one-request exporter,
and pins the separately audited post-buffer-marked CMS binary. It writes a
sealed capsule and build receipts; no target, exporter request, SAT query or
one-use claim is performed. Use `--validation-only` for the first build.
`icprog sat-target-publish` then archives all registered files, and
`icprog sat-target-replay` checks that archive as data without executing any
archived binary. Both require the capsule's external registration seal.

The first [validation result](result-v1/RESULT.md) passed complete archive
replay. The exporter built by that capsule also passed the three disclosed
source-byte and transport controls and a separate data-only audit. The
validation registration remains validation-only. The scientific freeze now
requires an independently replayed, exact-source validation publication and
three-control exporter audit. It rejects a different source commit, source
manifest or prepared-exporter binary. The one-use controller checks the
scientific publication and exact CMS source descriptor before it consumes the
claim; its frozen audit rechecks runtime receipts and the five online phases
as data. These new paths still need their own clean-commit validation capsule
and scientific publication. The other three source publications must also
exist before a fresh point card is drawn. No speedup follows from this build.
The [scientific sequence](SCIENTIFIC.md) names the exact admission checks and
one-use commands.
