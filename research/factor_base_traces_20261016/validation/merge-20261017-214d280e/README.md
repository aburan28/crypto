# Integration validation records

The receipt is `../../MERGE_VALIDATION_20261017_214d280e.md`.
Raw records are preserved once as `runs.tar.xz`. From this directory:

```sh
sha256sum -c runs.tar.xz.sha256
tar -xJf runs.tar.xz
```

The extracted `runs/` contains commands, original outcomes, full native
logs, source/input manifests and a per-member SHA-256 manifest. Names
prefixed `b6-before-merge-` belong only to the earlier b6 validation, not
the final source. Manifest paths denote extracted archive members.
Compiled binaries, credentials and chat transcripts are absent.

Deterministic repacking uses fixed metadata and the maintenance start date:

```sh
tar --sort=name --mtime='2026-10-17 00:00:00Z' --owner=0 --group=0 \
    --numeric-owner --mode='u+rwX,go+rX,go-w' --format=gnu -cf - runs \
  | xz -6 -T1 > runs-repacked.tar.xz
cmp runs.tar.xz runs-repacked.tar.xz
```

Previous committed run trees and receipts stay in their original form.
