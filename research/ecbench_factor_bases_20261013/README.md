# ecbench factor-base suite (2026-10-13)

Worked example of [`docs/ecbench/FACTOR_BASE_PROTOCOL.md`](../../docs/ecbench/FACTOR_BASE_PROTOCOL.md):
build dumps with `ecbench fb`, publish to S3 + Git, run `ic.pipeline` with the
same plug-in strings, check `fb_id` agreement.

## Catalog

See [`catalog.json`](catalog.json) (`ecbench.factor_base_catalog/v1`). Ten bases:

| id | plug-in | curve | store | git |
|:--|:--|:--|:--|:--|
| k17-orbit | `koblitz-orbit:divisor=0;1` | Koblitz a=1 n=17 | yes | dumps/ |
| k17-subspace | `binary-subspace:dimension=4` | same | yes | dumps/ |
| k17-compact | `compact-orbit-scan:columns=8,raw_x_cap=5000` | same | yes | dumps/ |
| k17-std | `koblitz-standard-subspace:dimension=4` | same | yes | dumps/ |
| k23-orbit | `koblitz-orbit:divisor=0;1` | Koblitz a=0 n=23 | yes | store only (411 KiB) |
| k23-compact | `compact-orbit-scan:columns=10,raw_x_cap=10000` | same | yes | dumps/ |
| k23-cyc-w1 | `cyclotomic:w=1` | Koblitz a=1 n=23 | yes | dumps/ |
| k11-cyc-w1 | `cyclotomic:w=1` | Koblitz a=1 n=11 | yes | dumps/ |
| p20-abscissa | `prime-abscissa:size=32` | 20-bit prime | yes | dumps/ |
| p20-glv | `glv-orbit:size=16` | same | yes | dumps/ |

Remote: `s3://crypto-autoresearcher/store` via `store/index.json`.

## Access

```bash
./target/release/ecbench store ls factor-bases/
./target/release/ecbench store pull store:factor-bases/icv1-f2m17-tm101-00378d4e/FB1hdd4118e0e739.json
./target/release/ecbench store verify --remote
```

## Run

```bash
cargo build --release --bin ecbench
./research/ecbench_factor_bases_20261013/run.sh
```

Or one session:

```bash
./target/release/ecbench run --wait --allow-busy \
  --spec research/ecbench_factor_bases_20261013/specs/koblitz17.json \
  --out research/ecbench_factor_bases_20261013/sessions/koblitz17
./target/release/ecbench verify --dir research/ecbench_factor_bases_20261013/sessions/koblitz17 --replay 8
./target/release/ecbench table --dir research/ecbench_factor_bases_20261013/sessions/koblitz17
```

## fb_id check

After a run, every IC arm's sealed `factor_base.fb_id` must equal the catalog
entry for that plug-in on that curve (see `check_fb_ids.py`).

## Claim hygiene

Public synthetic / known-answer curves only. Sessions are L0 operation-count
smokes on macOS unless re-run under L2 isolation. Not a performance or
sub-rho claim.
