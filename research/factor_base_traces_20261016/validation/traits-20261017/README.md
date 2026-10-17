# Native catalog repair validation

The receipt is `../../TRAITS_CI_REPAIR_20261017.md`. This maintenance does
not regenerate historical factor-base points or run a research campaign.

Raw records are preserved once in `runs.tar.xz`, with a SHA-256 sidecar.
From this directory:

```sh
sha256sum -c runs.tar.xz.sha256
tar -xJf runs.tar.xz
```

The extracted `runs/` contains the original CI failure, exact commands,
native logs, the original catalog fixture, source/input hashes, semantic
diff and per-member SHA-256 manifest. Paths in that manifest are archive
members. The final catalog is committed at `docs/curves/traits.json`;
the main-parent records are bound by the pinned Git commit in the receipt.
Compiled binaries and credentials are not included.

Deterministic repacking uses fixed metadata and the maintenance start date:

```sh
tar --sort=name --mtime='2026-10-17 00:00:00Z' --owner=0 --group=0 \
    --numeric-owner --mode='u+rwX,go+rX,go-w' --format=gnu -cf - runs \
  | xz -6 -T1 > runs-repacked.tar.xz
cmp runs.tar.xz runs-repacked.tar.xz
```

Previous committed evidence and validation receipts remain in their
original formats and locations. Nothing historical was repacked.
