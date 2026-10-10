# Sparse checkout for the harness

A full checkout is about 4.9 GB in 61,000 files; 4.5 GB of it is `research/`,
and 4.2 GB of that sits in 18 directories (tournament artifacts, `notes/`,
boolean-construction sweeps) that building ecbench, running `cargo test`, or
measuring a new session never reads. `scripts/sparse-checkout.sh` leaves them
off disk:

```sh
scripts/sparse-checkout.sh apply                     # the harness profile
scripts/sparse-checkout.sh status                    # what is on disk
scripts/sparse-checkout.sh add --path research/X/    # materialize more (end a directory with /)
scripts/sparse-checkout.sh patterns                  # print the rules, change nothing
scripts/sparse-checkout.sh disable                   # full checkout again
```

The profile is everything except those directories, plus every `research/`
path the Rust code needs, recomputed from the sources on each `apply`:
`include_str!`/`include_bytes!` targets anywhere (they must exist to compile)
and `research/` string literals in `src/` and `tests/` (read by `cargo test`).
On 2026-10-09 that was 42,413 of 60,752 files and 1.0 GB, every one of the 76
`include_*` targets present. New `research/<topic>_<date>/` directories, where
ecbench sessions land, are included by default, so `git add` accepts them.

## Downloading less

A sparse checkout of an existing clone frees disk only. For less network,
clone blobless and let the profile choose which blobs arrive:

```sh
scripts/sparse-checkout.sh clone https://github.com/aburan28/crypto DIR [--branch B]
```

In Claude Code on the web, `CRYPTO_SPARSE_PROFILE=harness` in the
environment's variables makes `.claude/hooks/session-start.sh` apply the
profile at session start.

## Off disk is not absent

An excluded file stays in `HEAD` and in the index; `git ls-files` lists it.
Before replaying, comparing against, or citing an existing session or note
under an excluded directory, `add --path` it. A path off disk is never
evidence that a record does not exist.

Examples that enumerate whole archives (`examples/n37_shared_rank_*` read
`research/notes/ecc2k130/`) need `add --path research/notes/ecc2k130/` before
they run. CI jobs keep their own per-job `sparse-checkout:` lists and are
unaffected.
