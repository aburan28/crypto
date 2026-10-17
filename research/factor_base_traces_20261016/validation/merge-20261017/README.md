# Integration validation records

The parent receipt is `../../INTEGRATION_VALIDATION_20261017.md`.
`check_preservation.py` is a read-only metadata/blob audit of the pinned
parents; it performs no curve arithmetic or measurement.

New raw validation records are packaged once as `runs.tar.xz`, following
AGENTS.md §13a. Verify and extract them from this directory:

```sh
sha256sum -c runs.tar.xz.sha256
tar -xJf runs.tar.xz
```

The extracted `runs/` contains commands, setup failures, native test/check
logs, source/input manifests, toolchain information and a per-member SHA-256
manifest. Manifest paths denote extracted archive members, not loose Git
paths. No credentials or chat transcript are included.

Repacking the extracted tree uses the validation start date and fixed
metadata, so it yields identical archive bytes:

```sh
tar --sort=name --mtime='2026-10-17 00:00:00Z' --owner=0 --group=0 \
    --numeric-owner --mode='u+rwX,go+rX,go-w' --format=gnu -cf - runs \
  | xz -6 -T1 > runs-repacked.tar.xz
cmp runs.tar.xz runs-repacked.tar.xz
```

The earlier `catalog-20261016/` and `merge-20261016/` records stay in their
committed formats; no previous evidence tree was repacked or rewritten.
