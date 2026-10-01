# Four-window screen source lock

This change implements the [preregistered screen](PROTOCOL.md) on the registered
`icv1-f2m41-tm2308219-7f48b14a` rung. It contains no new public targets,
scored window-1/2/3 headers, charged timing, or conclusion about a faster
discrete-log method. The frozen input PR and one-shot outcome PR are separate
dependencies. The batch screen is a **pre-selector opportunity diagnostic**;
it does not replace the repository's one-target IC-versus-rho acceptance gate.

`build_frozen.py` materializes the 20 pinned v2 source files, copies in the
merged window generator and frozen Cargo.lock, checks every source hash, and
builds the generator, compact solver, and normal-basis rho in one synthetic
checkout. It writes a materialization record, compiler versions, build logs,
binary hashes, and `BUILD_RECEIPT.json`. The source-only CI build uses the dev
profile; its generator must reproduce the frozen toy and archived n=41 window-0
headers byte for byte. That control is not a scored run.

After this source-lock PR merges, use its **merge commit** as the immutable
source anchor. In a full checkout with all previously tracked point corpora
present, create a new input-freeze branch and run:

```sh
python3 research/notes/ecc2k130/base_window_screen_20261001/prepare_inputs.py \
  --source-merge-commit <40-hex-merge-commit>
python3 research/notes/ecc2k130/base_window_screen_20261001/verify_inputs.py \
  --out research/notes/ecc2k130/base_window_screen_20261001/INPUT_RECEIPT.json
```

Review and merge the five point-only files, separate known-scalar fixtures,
`FROZEN.json`, and independent replay receipt before building or running a
scored arm. The preparer checks that every source file matches the merge
commit, inventories all earlier tracked point corpora from the preparation
HEAD's immutable Git tree, and deterministically
derives orbit-disjoint Q. The independent verifier replays each accepted and
rejected candidate and all 5,120 group equations. Neither command samples a
base-window or method cost.

On the Linux x86-64 timing host, after the input freeze merges, use an empty
build directory and a single reserved physical CPU plus its SMT sibling:

```sh
python3 research/notes/ecc2k130/base_window_screen_20261001/build_frozen.py \
  --out <absolute-build-directory> --profile release --offline
python3 tools/isolated_bench.py reserve \
  --cpus <cpu,smt-sibling> --settle 2 --period 1 \
  --max-other-cpu 0.10 --label base-window-screen-n41 \
  --out <absolute-run-directory>/isolation.jsonl -- \
  python3 research/notes/ecc2k130/base_window_screen_20261001/run_screen.py \
    --source-root <absolute-build-directory>/source \
    --generator <absolute-build-directory>/target/release/examples/koblitz_base_window \
    --compact <absolute-build-directory>/target/release/examples/koblitz_orbit_dlp_s3_batch \
    --rho <absolute-build-directory>/target/release/examples/koblitz_rho_batch_ks_v3 \
    --materialization <absolute-build-directory>/materialization.json \
    --build-receipt <absolute-build-directory>/BUILD_RECEIPT.json \
    --out <absolute-run-directory> --cpu <cpu>
```

The runner creates the run directory and never opens the known-scalar files.
It charges a fresh generator plus a fresh compact process for every window
arm, and one rho process on the same Q in each block. Keep the run directory,
including partial reports after a timeout or failure, without overwriting it.
Save a preflight refusal separately; do not convert it into a successful or
negative timing outcome. If crates are unavailable for an offline release
build, resolve that **before** entering the measured reserve window and record
the dependency acquisition in the build receipt.

Replay and classify the immutable run with:

```sh
python3 research/notes/ecc2k130/base_window_screen_20261001/verify_screen.py \
  --run-dir <absolute-run-directory> --out <absolute-replay.json>
python3 research/notes/ecc2k130/base_window_screen_20261001/analyze_screen.py \
  --run-dir <absolute-run-directory> --replay <absolute-replay.json> \
  --out <absolute-analysis.json>
```

The independent replay checks all frozen source/input hashes, the generator
scan and all orbit members, every rank transition and four-point witness,
and every recovered compact and rho scalar. The analyzer then charges child
user-plus-system CPU and applies the frozen A/A, isolation, fixed-window, and
hindsight rules. A passing dev build or source parity does not supply a timing
result. A favorable hindsight lower envelope is not a selection policy or
an ECC2K-130 speedup. Publish the raw run, replay, and decision in a focused
outcome PR; update the canonical scoreboard only after independent review.
