# Native optimization fixtures and later output verification

The generator froze **three training points, then five holdout points** with
SplitMix64 seed `202610128301`. All eight are nonidentity points on the exact
legacy field/model, belong to the prime subgroup and pass native `[d]G=Q`
replay. Their signed Frobenius orbits are distinct from the prior fixed Q,
the retained K600 factor base and one another. The independent fixture replay
also reproduced the deterministic draw/rejection stream and file hashes.

The exact model is `icv1-f2m83-t6151469093347-9eba6cd5`, with field low terms
`[0,2,4,7]`, subgroup prime `8569786107849059` and the retained generator.
The historical `cdcc5432` alias belongs to another modulus and is not used
for these fixtures.

Canonical frozen artifacts are under [../fixtures/](../fixtures/), with
byte-for-byte provenance copies in this directory's `fixtures/` subdirectory:

- `training_1.jsonl` through `training_3.jsonl`: one point per file.
- `holdout_1.jsonl` through `holdout_5.jsonl`: one point per file.
- `holdout_5_targets.jsonl`: the five point-only holdout rows in order.
- `validation_sidecar.json`: separate validation-only scalar values.
- `generator_receipt.json`: seed, law, counts, exact identity, points and hashes.
- `artifact_hashes.sha256`: hashes of the generated artifacts.

Measured children receive only point files. Never pass the scalar sidecar to
a measured solver. The proposal scalar law is uniform on `1..r-1` using
64-bit rejection sampling. Accepted points are conditioned on the specified
orbit exclusions and selected without repeated signed orbits; they are not
claimed to be an unconditioned IID target sample. The eight accepted points
used eight draws, with zero range or orbit rejections in this seed.

Generation receipt SHA-256:

```text
c93edc32a77411443229191e2db4cc3eeb7e99976b7f79a8a753ba103bce9909
```

`fixture_replay_receipt.json` reports `PASS_FIXTURES_ONLY`, with all eight
fixtures validated and zero invalid rows. Its SHA-256 is:

```text
fb263efde55e8177466ced20df6f09bd3b74c60b7c7505b62f91d6a650860315
```

Native sources are `generate_targets.cpp`, `check_optimization.cpp`,
`native_helpers.hpp` and the earlier standalone arithmetic `check.cpp`.
The helper includes that arithmetic source after renaming its archived
entry point; it never invokes the old archived checker/report code and
imports no Rust solver. Compiled executables remain outside the repository;
their original hashes are in `executable_hashes.sha256`.

```sh
g++ -std=c++17 -O2 -Wall -Wextra -Werror generate_targets.cpp -o /tmp/n83-generate
g++ -std=c++17 -O2 -Wall -Wextra -Werror check_optimization.cpp -o /tmp/n83-verify
/tmp/n83-generate ARCHIVED_FIXTURE_BASE_DIR NEW_FIXTURE_DIR
/tmp/n83-verify --fixtures-only ARCHIVED_FIXTURE_BASE_DIR NEW_FIXTURE_DIR NEW_RECEIPT
/tmp/n83-verify ARCHIVED_FIXTURE_BASE_DIR NEW_FIXTURE_DIR COMPLETED_PANEL_DIR NEW_RECEIPT
```

`ARCHIVED_FIXTURE_BASE_DIR` is
`experiments/koblitz-single-target-n83-20261006`. Generation refuses an
existing output directory and creates files using `O_EXCL`. Both checker
modes refuse an existing receipt and use exclusive creation. The fixture
files and original records must not be overwritten when replaying; use a
fresh directory/receipt path.

Run output verification only after all timed children finish. It recursively
reads strong `producer.json` and IC `target.jsonl`/`targets.jsonl`, checks
every scalar answer against its exact frozen public point, verifies each
IC four-point witness against the independently checked 99,600-point base,
and preserves cap/failure/missing-answer process dispositions. The scalar
sidecar is compared after the recovered answer's point replay.

Output verification is compiled and ready; no optimization benchmark or its
completed output replay is claimed here. Full factor-log/rank transcript,
performance reproducibility and second-host/isolation gates remain separate
obligations. The parent task owns the protocol, runner, timed records,
statistics and final study report.

The fixture-generation [report](GENERATION_REPORT.md) and editable
[diagram](generation.svg) describe this freeze without a timing claim.
Conductor tooling was unavailable; work stayed within the explicitly
delegated optimization/verification and optimization/fixtures paths. No
old evidence was edited and no git commit or external write was made here.
