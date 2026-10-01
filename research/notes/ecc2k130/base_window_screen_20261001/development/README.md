# Base-window generator source controls

These are **source-only** controls. No held-out Q, candidate windows 1–3,
rank/target measurements, or rho processes are involved. The committed
`SOURCE_CONTROL.json` is reproduced by:

```sh
cargo build --locked --example koblitz_base_window
python3 research/notes/ecc2k130/base_window_screen_20261001/test_generator.py \
  --generator target/debug/examples/koblitz_base_window \
  --out /tmp/base-window-source-control.json
cmp /tmp/base-window-source-control.json \
  research/notes/ecc2k130/base_window_screen_20261001/development/SOURCE_CONTROL.json
```

The toy reference header was emitted with the **frozen** compact source
`702a0a05709bc14bc10bafbb08edb2a2f5f86794e970d84c317da3bc65cdcf38`,
materialized from the v2 freeze at
`2c207b41a3be3552c1754ee78540f13be258a105`. Its source tree used the
frozen lockfile SHA-256
`7d671f48c2da93f133d98802e80f858d1d9ea3b86996f7037f758990e1566627`.
The input was a one-line legacy scalar fixture `1\n`, and the command was:

```sh
KIC_DUMP_BASE=/tmp/n13_k3_frozen_construct.base.jsonl \
  target/debug/examples/koblitz_orbit_dlp_s3_batch \
  construct:13:0:3 /tmp/toy-scalar-1.txt 7 /tmp/toy-target.jsonl
```

The frozen process completed rank 3 and solved its one toy target. The
committed header SHA-256 is
`69fa0b23ab0604ecc7259dbe88770bf0514bed26168cb16eaf9b86d23fa41545`;
the new generator reproduces its bytes exactly. The n41 control compares
against the committed v2 cold-panel archive
`evidence_run_36803331080/raw/n41_L1024.tar.gz`, SHA-256
`261b5f3e7b209c2163d807dc78554b96b80873eda052405156b1178b29dce4a4`.
The independent verifier replays every scan decision and every selected
signed-Frobenius orbit member with separate Python field and curve arithmetic.
The raw-x cap failure, existing-output refusal, and a deliberately corrupted
scan are also tested.
