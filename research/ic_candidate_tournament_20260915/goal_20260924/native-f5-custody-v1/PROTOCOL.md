# Full native F5 capsule custody

The native F5 controller builds from a frozen offline snapshot, but its compact
build-validation bundle deliberately omits most source/dependency/binary bytes.
It cannot establish custody of an actual executable registration. This follow-on
adds native publication and portable data replay of the full capsule before its
sole scientific dispatch. It registers no job, generates no target and launches
no solver. The disclosed synthetic n17 scope and prior closed registrations
remain unchanged.

`icprog f5-control-publish-custody` accepts only an unconsumed executable-mode
capsule, checks its complete live inventory and mathematical-only input, copies
its six exact sidecars with create-only writes, and audits a separately packed
plain USTAR gzip archive. It rechecks the capsule and absence of a consumed claim
after copying. The full archive contains only `immutable/`, `config.json`,
`preparation.json`, `job.json`, `host-context.json`, `registration.json` and
`seal.json`. Build caches and later execution claims are excluded.

`icprog f5-control-replay-custody` requires the externally published registration
SHA-256. It checks canonical registration identity, every sidecar, prepared
geometry/logs and mathematical-only worker input, source and binary inventory,
both executable digests, compressed archive bytes and every decompressed member.
Compressed input is bounded at 96 MiB; decompressed data is bounded at 320 MiB.
It reads archive bytes as data and extracts or executes nothing. A custody pass
establishes retained bytes only: execution admission, native runtime qualification,
fresh-target qualification, headline eligibility and promotion remain false;
online time and speedup remain unknown.

The archive parser is shared with the existing SAT registration auditor. Its
original checksum/USTAR/path/member/trailer rejection controls remain applicable.
No historical capsule or registration is rewritten. The original frozen SAT
checker and all consumed invocations stay closed. Portable postexecution data
tools do not replace either family's original preexecution-frozen auditor.

Use the native busy wrapper. Paths and the seal below are placeholders, not an
actual executable registration or authorization to dispatch a restored capsule:

```sh
# Create a new empty publication directory first. Thin shell packs data only.
COPYFILE_DISABLE=1 tar --format ustar -czf NEW_PUBLICATION/capsule.tar.gz \
  -C NEW_CAPSULE immutable config.json preparation.json job.json \
  host-context.json registration.json seal.json
/tmp/ic-native-busy busy -- target/debug/icprog f5-control-publish-custody \
  --capsule NEW_CAPSULE --publication NEW_PUBLICATION --out NEW_PUBLICATION_RECEIPT
/tmp/ic-native-busy busy -- target/debug/icprog f5-control-replay-custody \
  --publication NEW_PUBLICATION --registration-sha256 PUBLISHED_SEAL \
  --out NEW_REPLAY_RECEIPT
```

Publish the actual archive and sidecars in a focused preregistration PR. Require
reviewed source, replayable durable publication and all applicable exact-head CI
before the sole claim-consuming invocation. A local path alone is not publication.
Original audit then uses the frozen checker and the same external seal. Preserve
success, incompleteness, timeout, crash and audit failure without retry or seed
replacement. The full goal still requires new native natural-yield evidence and
a new frozen one-target comparison with the incumbent and matched rho.

Tests use placeholder archive bytes and never run a worker. They check positive
custody separately from execution admission, external-seal mismatch, create-only
output, altered/resealed sidecars, omitted evidence, truncated archives and
prohibited execution claims. The existing SAT publication is replayed after the
shared-parser refactor. This is a reproducibility change, not an IC speed result.
