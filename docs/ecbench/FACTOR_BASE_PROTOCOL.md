# Factor-base protocol for ecbench

One way to name, build, store, and run a factor base so a Git checkout, a
CI job, and an S3-backed worker all resolve the same bytes and the same
`FB1` identity. This protocol prices nothing; it only makes the inventory
and the `ic.pipeline` plug-in string agree.

## Roles of each place

| place | holds | role |
|:--|:--|:--|
| **Plug-in spec** | e.g. `koblitz-orbit:divisor=0;1` | what `ic.pipeline` rebuilds at run time; the canonical constructor |
| **Dump** (`ecbench.factor_base_dump/v1`) | curve facts, `FactorBaseFacts` (`fb_id`, point digests), every point | inventory / database / claim identity; not the live matrix input today |
| **Git** | dump JSON when small; always the store index | reviewable pointer + small fixtures |
| **S3 store** (`s3://crypto-autoresearcher/store`) | dump bytes (any size) | durable content-addressed blob |
| **`store/index.json`** | name → sha256, bytes, provenance | Git-tracked manifest (`docs/store/README.md`) |

`ic.pipeline` still builds from the plug-in string. The dump's `fb_id` must
match the id the run records for that string on that curve; that is how a
stored dump certifies a session without becoming a second constructor.

## Naming

Store names (and dump basenames) follow `docs/store/README.md`:

```
factor-bases/<icv1-slug>/<FB1-id>.json
```

Example:

```
factor-bases/icv1-f2m17-tm101-00378d4e/FB1hdd4118e0e739.json
```

Do not invent parallel slug spellings. Curve constructors must expand to a
registered ICV1 slug (`ecbench plan` prints `registered=true`).

## Build

```bash
cargo build --release --bin ecbench

./target/release/ecbench fb \
  --curve '{"kind":"koblitz","a":1,"n":17}' \
  --factor-base 'koblitz-orbit:divisor=0;1' \
  --out /tmp/fb.json
```

The tool prints `FB1h…` and the path. The JSON schema is
`ecbench.factor_base_dump/v1`.

Supported plug-in families (binary / Koblitz unless noted):

| family | example | notes |
|:--|:--|:--|
| `koblitz-orbit` | `koblitz-orbit:divisor=0;1` | signed Frobenius fold |
| `binary-subspace` | `binary-subspace:dimension=4` | abscissa subspace |
| `compact-orbit-scan` | `compact-orbit-scan:columns=10,raw_x_cap=10000` | raw-x scan + fold |
| `cyclotomic` | `cyclotomic:w=1` | needs `2n+1` prime, `ord(2)=n` (e.g. n=11,23) |
| `onb-weight` | `onb-weight:w=1[,columns=N]` | ONB Hamming-weight-`w` abscissae; same fold as cyclotomic; odd `n<64` |
| `koblitz-standard-subspace` | `koblitz-standard-subspace:dimension=4` | PDP / standard basis (binary plug-in) |
| `prime-abscissa` | `prime-abscissa:size=32` | prime curves |
| `glv-orbit` | `glv-orbit:size=16` | prime curves with automorphism fold |
| `prime-kummer` | `prime-kummer:t=9` | prime curves, `2^t \| p − 1`; PKM's Kummer tower |
| `prime-isogeny` | `prime-isogeny:t=9,seed=0` | prime curves, any `p`; PKM's isogeny tower, auxiliary-curve search charged to the factor-base phase |

## Publish (Git + S3)

**Small dumps (≤ 256 KiB)** also live in the research tree next to the
catalog for offline reading. **Every dump** is pushed to the store so the
index is the durable pointer:

```bash
FB_ID=$(jq -r .factor_base.fb_id /tmp/fb.json)
SLUG=$(jq -r .curve.slug /tmp/fb.json)
NAME="factor-bases/${SLUG}/${FB_ID}.json"

./target/release/ecbench store push /tmp/fb.json \
  --name "$NAME" \
  --description "ecbench dump of koblitz-orbit:divisor=0;1 on ${SLUG}" \
  --provenance "ecbench fb --curve '…' --factor-base 'koblitz-orbit:divisor=0;1' (commit $(git rev-parse --short HEAD))"

git add store/index.json
# and, when ≤ 256 KiB:
# cp /tmp/fb.json research/ecbench_factor_bases_YYYYMMDD/dumps/${FB_ID}.json
```

Access anywhere:

```bash
./target/release/ecbench store pull "store:${NAME}"          # prints cache path
./target/release/ecbench store info "store:${NAME}"
./target/release/ecbench store verify --remote "$NAME"
```

`crypto store …` is the same command. A plain path argument that accepts
`store:NAME` (today: `ic … --factor-base`) resolves through
`crypto_lib::store::resolve_path`.

## Load into the database

```bash
./target/release/ecbench db sql path/to/dump.json | sqlite3 -bail ecbench.db
```

Sessions that used the same plug-in on the same curve join on `fb_id`.

## Run ecbench with that base

Write an `ecbench.spec/v1` whose IC arm repeats the **same** plug-in string
used to build the dump. Example arm:

```json
{"name": "ic-orbit", "role": "candidate", "method": {"id": "ic.pipeline", "params": {
  "factor_base": "koblitz-orbit:divisor=0;1",
  "oracle": "mitm-frobenius:m=2"
}}}
```

```bash
./target/release/ecbench plan --spec SPEC.json
./target/release/ecbench run --wait --allow-busy --spec SPEC.json --out OUT_DIR
./target/release/ecbench verify --dir OUT_DIR --replay 8 --out OUT_DIR.audit.json
./target/release/ecbench table --dir OUT_DIR
```

`mitm-frobenius` needs a Frobenius-folded base (`koblitz-orbit`,
`compact-orbit-scan`, `cyclotomic`, …). Plain `binary-subspace` and
`koblitz-standard-subspace` use `mitm:m=2` (or another non-folding oracle);
pairing them with `mitm-frobenius` fails with
`this instance has no Koblitz structure to fold by`.

After the run, check the sealed record's factor-base id equals the dump:

```bash
jq -r '.factor_base.fb_id' dumps/FB1h….json
jq -r 'select(.arm=="ic-orbit") | .factor_base.fb_id' OUT_DIR/records.jsonl | sort -u
```

A mismatch means the constructor or curve slug drifted; do not cite the
session against that dump.

## FBM1 — measurement identity

FB1 names a point set. **FBM1** names the measurement configuration that
produced yield / DoR / cost numbers against that set:

```
curve (ICV1 slug) × FB1 × oracle × solver × linalg × instrument
```

Results (S, yield, solving degree, Sage D\*, …) attach to an FBM1 id; they
are never part of the hash (same rule as ICV1). Schema:
`ecbench.factor_base_measurement/v1`. Short id: `FBM1h` + 12 hex.

```bash
./target/release/ecbench fbm \
  --curve '{"kind":"koblitz","a":1,"n":17}' \
  --factor-base 'binary-subspace:dimension=4' \
  --oracle 'descent-algebraic:m=2' \
  --solver 'f4-f2' \
  --out /tmp/fbm.json
```

IC records written by `ecbench run` carry `fbm_id` next to `factor_base`
when the run used a factor base. The database tables
`factor_base_measurements` / `factor_base_measurement_cells` and the view
`fb_measurement_ledger` index them (`docs/ecbench/schema.sql` v3).

**Cite yield or DoR only with both `fb_id` and `fbm_id`.** Orphan Sage DoR
JSON (no matching FB1 dump / catalog entry) is refused by
`research/ecbench_factor_bases_20261013/join_sage_dor.py`.

Instrument tags:

| tag | meaning |
|:--|:--|
| `ecbench.ic_pipeline/v1` | native `ic.pipeline` session |
| `sage.factor_base_dor/v1` | Sage `factor_base_dor.sage` cell joined by FB1 |

## Algebraic panel

Table oracles (`mitm`, `mitm-frobenius`) fill yield and S. Algebraic
oracles (`descent-algebraic` / counted) fill `run_solver` (DoR proxies).
Keep mitm arms as yield/cost controls next to algebraic arms on the same
FB1. Specs live under `research/ecbench_factor_bases_20261013/specs/`
(`*-algebraic.json`, `*-onb.json`).

## Catalog

A research round that freezes several bases keeps a catalog JSON (schema
`ecbench.factor_base_catalog/v1`) listing for each entry: plug-in string,
curve constructor, `fb_id`, store name, optional Git relative path, and
the oracle used in the companion ecbench spec. Grow the catalog with
every new family. CI runs `check_fb_ids.py` (and `check_fbm_ids.py` when
sessions carry `fbm_id`). The worked example is
`research/ecbench_factor_bases_20261013/`.

## Git vs S3 decision

| size | Git | Store |
|---:|:--|:--|
| ≤ 256 KiB | dump under `research/…/dumps/` | always (`store push`) |
| > 256 KiB | catalog + index entry only | always |

Never commit a dump without either (a) the bytes in Git or (b) a
`store/index.json` entry whose digest verifies. Prefer both for smoke
sizes.

## Related

- Store mechanics: [`docs/store/README.md`](../store/README.md)
- Harness: [`docs/ecbench/README.md`](README.md)
- Claim / FB1 identity: README §4 and `ecbench claim`
