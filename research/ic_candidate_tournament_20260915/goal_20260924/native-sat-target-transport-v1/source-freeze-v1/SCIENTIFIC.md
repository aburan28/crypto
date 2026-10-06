# SAT source publication sequence (no point card yet)

The first validation-only result in [`result-v1`](result-v1/RESULT.md)
predates the scientific controller. It demonstrates the build and parity
paths, but its source commit cannot authorize the newer controller. Start a
new validation capsule from one clean committed source snapshot that already
contains `sat-target-execute` and `sat-target-audit`.

1. Freeze a new target-free `sat-target-freeze --validation-only` capsule with
   this directory's `config.json`, `preparation-binding.json`, and
   `host-context.json`, the exact marked-CMS binary, and explicit `cargo` and
   `rustc` paths. The command runs offline build tools and zero target queries.
   Save its external registration seal. `sat-target-publish` and
   `sat-target-replay` with that seal must agree on every archived file.
2. Run all three fixed `ic_exporter_prestart_probe` controls against the
   prepared exporter **from that capsule** and the accepted exporter. Preserve
   each raw control directory. Use `sat-exporter-audit` with the validation
   publication, seal, raw controls, exact probe binary, and a new output path.
   The audit must report `PASS_DISCLOSED_EXACT_BUILT_EXPORTER_PARITY` and the
   same source manifest/archive/exporter hashes as the validation publication.
3. Without changing or recommitting source, run a separate
   `sat-target-freeze` without `--validation-only`. Supply all five extra
   arguments: `--validation-publication`,
   `--validation-registration-sha256`, `--exporter-controls`,
   `--exporter-probe`, and `--exporter-audit`. Use
   `host-context-scientific.json`. This freeze rejects a changed Git commit,
   source inventory, parity result, or prepared-exporter binary. Save its new
   external seal; publish and data-only replay it with `--scientific` before
   any public point exists.
4. The campaign must independently inspect this full archive and the other
   three source publications, then publish all four strict source descriptors.
   Only then may it draw one scalar-free `hash_to_subgroup_v1` point card. The
   card's `cms` source hash must be the exact descriptor byte hash. Its point,
   resource envelope, and arm order remain fixed for all four one-use runs.
5. Execute only the `icprog` binary inside the new scientific capsule with
   `sat-target-execute --capsule ... --publication ...
   --source-descriptor ... --card ... --execution NEW_DIR
   --registration-sha256 ...`. It consumes the registration before launching
   one worker. A failure is terminal; never retry or restore the capsule.
   Use the same frozen `icprog` for `sat-target-inspect` and `sat-target-audit`
   with outputs outside the capsule and execution. The audit replays original
   source/model files and the five online phases without a solver call.

Every validation, publication and disclosed parity result remains a stage
control. Even a verified one-target recovery on this ordinary host would be
exploratory; a CPU speedup needs a same-point rho solve and the host-level
isolation receipt required by `AGENTS.md`.
