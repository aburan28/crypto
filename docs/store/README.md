# The artifact store: large files without Git

Factor bases, relation tables, DP archives and raw evidence that are too big
for Git live in object storage and are fetched when something needs them.
Code names a file, not a place:

```
ic run --factor-base store:factor-bases/icv1-f2m41-tm2308219-7f48b14a/fb1-0123abcd.json
```

The store turns that name into a verified local path, downloading it once
into a content-addressed cache. A small index tracked in Git records what
each name means, so a checkout, a CI job, a worker pod and the
`crypto-autoresearcher` repository all resolve the same name to the same
bytes.

The store moves bytes and checks them. It prices nothing and makes no
performance claim; the evidence rules in `AGENTS.md` apply to what is built
on the bytes.

## Quick start

```bash
# once per machine: the AWS CLI with credentials that can use the bucket
aws sts get-caller-identity

# publish a file (uploads, caches, records it in store/index.json)
crypto store push fb.bin \
  --name factor-bases/icv1-f2m41-tm2308219-7f48b14a/fb1-0123abcd.bin \
  --description "subspace d=8 factor base" \
  --provenance "ic search --curve ... --out fb.bin (commit abc1234)"
git add store/index.json && git commit      # the index is the manifest

# anywhere else: fetch on demand
crypto store pull store:factor-bases/icv1-f2m41-.../fb1-0123abcd.bin   # prints the cached path
crypto store pull <name> --out ./fb.bin                                # also copy it out
```

`ecbench store …` is the same command. A path argument that accepts a
`store:NAME` or `sha256:DIGEST` reference needs no other change; today that
is `ic … --factor-base`, and any code can adopt it with one call:

```rust
let path = crypto_lib::store::resolve_path(&arg)?;   // a plain path is returned unchanged
```

## How it works

| piece | what it is |
|:--|:--|
| **object** | immutable bytes, keyed by their SHA-256: `sha256/ab/abcdef…` |
| **index** | `store/index.json`, tracked in Git: name → sha256, bytes, description, tags, provenance, superseded digests |
| **cache** | `$CRYPTO_STORE_CACHE`, default `~/.cache/crypto-store`; verified objects, read-only, written atomically |
| **remote** | `s3://bucket/prefix` through the `aws` CLI, or `file:///directory` |

Why content addressing: a cached copy never goes stale, a remote object is
never overwritten with different bytes, and every download is hashed against
the digest the index promised, so a wrong or damaged object is refused and
never reaches the cache.

The index is the manifest `AGENTS.md` asks for when output is too large for
Git: content hash, byte count, durable location (the remote in
`store/config.json`), and the command that produced it (`--provenance`).
A name that already holds other content is not repointed without
`--replace`, and the old digest is kept in `superseded`, so a replaced
figure stays visible as the "was".

## Commands

| command | does |
|:--|:--|
| `store status` | the index, cache and remote this checkout resolves to, and whether the remote answers |
| `store ls [prefix]` | indexed names, with whether each is cached |
| `store info <spec>` | one entry in full |
| `store push <file> --name N` | upload (skipped if the remote already holds the bytes), cache, index |
| `store pull <spec> [--out P]` | make available locally, print the path |
| `store verify [names] [--remote]` | re-hash cached objects; with `--remote`, confirm each object exists with the indexed size and digest |
| `store gc [--apply]` | drop cached objects the index no longer names (dry run by default) |
| `store rm <name>` | remove a name from the index |

A `<spec>` is a name, `store:NAME`, or `sha256:DIGEST` (a full digest, or a
prefix of 8 or more digits that the index resolves).

The remote is **never deleted from**: `rm` edits the index, `gc` prunes the
local cache. Deleting remote objects is a decision for the bucket's owner.

## Configuration

Precedence is environment, then `store/config.json`, then the default.

| setting | environment | `store/config.json` | default |
|:--|:--|:--|:--|
| remote | `CRYPTO_STORE_REMOTE` | `remote` | none (the committed config names `s3://crypto-autoresearcher/store`) |
| endpoint (MinIO, LocalStack, …) | `CRYPTO_STORE_ENDPOINT` | `endpoint` | AWS |
| cache directory | `CRYPTO_STORE_CACHE` | | `$XDG_CACHE_HOME/crypto-store`, else `~/.cache/crypto-store` |
| index file | `CRYPTO_STORE_INDEX` | | `store/index.json` at the repository root |
| `aws` executable | `CRYPTO_STORE_AWS` | | `aws` |

A local directory works as a remote (`CRYPTO_STORE_REMOTE=file:///srv/store`),
which is how a lab without S3, a shared volume, or a test uses the store.

## S3 details

Transport is the AWS CLI, as in `ecc2k-merge` and `isogeny_walk::store`, so
credentials, regions, SSO and retries are the CLI's and this crate links no
HTTP stack. Uploads follow the `isogeny_walk::store` conventions:

* **Write-once.** `put-object --if-none-match '*'`: an existing key is never
  overwritten. If a concurrent writer stored the key first, the store accepts
  it only after checking the size and digest match.
* **Checked on arrival.** The upload declares `--checksum-sha256`, which S3
  verifies against the bytes it receives, and the reply's checksum is
  compared again.
* **Large objects.** Past 4 GiB a single PUT is not possible, so the upload
  is multipart through `aws s3 cp`, with the digest in user metadata. S3
  cannot check a whole-object SHA-256 across parts, so the guarantee for
  those objects rests on the hash taken on download.
* **Checked on download.** Every object is hashed as it enters the cache, and
  a mismatch is discarded.

The IAM permissions needed are `s3:GetObject`, `s3:PutObject` and
`s3:ListBucket` on the bucket and prefix. S3 answers `403` rather than `404`
for an absent key when `s3:ListBucket` is missing; the error says so.

## Naming

Names are path-like (`[A-Za-z0-9._/+=-]`, no `..`). Use a namespace that says
what the file is and which curve it belongs to, with the ICV1 slug (AGENTS.md
§11) for curves:

```
factor-bases/<icv1-slug>/<FB1-id>.<ext>
tables/<icv1-slug>/<description>.<ext>
evidence/<topic>_<date>/<file>
```

Do not rename frozen evidence to fit: the store name is a pointer to the
bytes, the committed file keeps the name it was written with.

## For evidence

When a result's raw output is too large for Git, push it and commit the
index entry with the PR. The manifest an evidence directory needs (hash,
byte count, location, extraction command) is then the entry's `sha256`,
`bytes`, the remote in `store/config.json`, and `provenance`. A reader
reproduces the bytes with `crypto store pull sha256:<digest>`.

## What this does not do

* **Not tested against live S3 from the session that wrote it.** The
  argument vectors are unit-tested, and the whole flow runs end to end
  through the real subprocess boundary against a stand-in `aws` that models
  404s, `If-None-Match` and checksum rejection (`tests/store_cli.rs`). The
  first push to the real bucket is the first live check; `store status`
  reports whether the remote answers.
* **No migration of existing large files.** Files already in Git history stay
  there; moving them out is a separate, history-visible decision.
* **No remote deletion, versioning policy or lifecycle rules.** Those belong
  to the bucket.
* **No database.** An RDS-backed catalogue of artifacts (and of run yields)
  is the next step; the index here is its source of truth until then.

## Related

* `src/cryptanalysis/isogeny_walk/store.rs`: the run-level publisher for walk
  outputs (a completion marker over many objects). It shares the write-once
  and checksum conventions; this store is for named files of any kind.
* `ecc2k130/aws/`: the ECC2K-130 campaign's own bucket and tooling.
* `taskq/`: workers need the `aws` CLI (or a `file://` remote on a shared
  volume) to resolve `store:` references.
