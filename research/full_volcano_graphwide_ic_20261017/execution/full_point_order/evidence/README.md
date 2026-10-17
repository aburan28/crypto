# Full selected-point-order C evidence

`runs.tar.xz` is the final deterministic delivery (9,017,364 bytes). It holds
all four raw cohorts, exact executable and tests, source, pinned dependency
digests, C prefix and local certificates, phase records, truth, isolation logs,
native summaries, and inspected reports. Replays remain separate populations.

```sh
sha256sum -c runs.tar.xz.sha256
mkdir runs
tar -xJf runs.tar.xz -C runs
(cd runs && sha256sum -c MANIFEST.sha256)
```

All manifest entries passed after extraction. The extracted measured executable
re-ran the original two-target control and matched every non-clock record,
fixture and truth field. Repacking the extracted tree gave the same outer
SHA-256; receipts are beside this README. This is same-implementation replay,
not independent arithmetic validation. Parent README and frozen protocol state
the smooth-order, unknown-conductor, timing and full-population limits.

Loose native summaries/CSV, isolation journals and replay receipts are compact
review pointers; raw records remain in the archive. `delivery.json` records
local hashes and Drive name/size/parent readbacks, not remote content hashing.
The B cache-eviction receipt concerns only a reproducible compressed duplicate;
its authoritative 9.3-GB raw graph remains local and its prior Drive parts were
not deleted. No S3 upload is claimed.
