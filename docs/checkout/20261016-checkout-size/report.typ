#set page(paper: "a4", margin: 2cm)
#set text(size: 10pt)
= Checkout size by recipe, 2026-10-16

`main` at `f25d5e520`, cloned once per recipe from github.com through a Claude Code cloud container's HTTPS proxy (shared network, no isolation). Recipes: `docs/sparse-checkout.md`. Values: `docs/checkout/20261016-checkout-size/data.json`.

#figure(image("checkout-size.svg", width: 100%), caption: [Bytes on disk after the clone (working tree and `.git`), files on disk and wall time; one clone per recipe at `f25d5e520`.])

#table(columns: 6, align: (left, right, right, right, right, left),
  [*recipe*], [*files*], [*tree (bytes)*], [*.git (bytes)*], [*time (s)*], [*`cargo check --all-targets`*],
  [git clone], [100,218], [6,651,988,041], [4,061,787,984], [190], [not run],
  [blobless, harness], [81,435], [2,521,984,680], [749,208,026], [46 after minimal], [pass],
  [blobless, minimal], [26,078], [458,328,475], [140,696,247], [33], [pass],
  [blobless, bare cone], [3,434], [220,299,783], [62,890,969], [9, 8], [fail; lib and bins pass],
)

All values measured (`du -sb`; wall seconds around the clone). One sample per recipe, two for the bare cone; times show order of magnitude only. The harness row applied `--profile harness` to the minimal clone.

Derived from `git ls-tree -r -l f25d5e520`: 829 of 100,218 tracked files are over 1 MiB and hold 4,823,311,258 of 6,651,988,041 bytes (72.5%), the part `blob_limit_kb: 1024` leaves out of a cloud session's clone.

The bare cone of `src docs scripts tools benches examples tests` fails `cargo check --all-targets` because four Cargo targets are declared in `experiments/`; the minimal profile adds those directories and every `include_str!`, `include_bytes!` and `#[path]` target, and passes (160 s cold, 30.8 s warm).
