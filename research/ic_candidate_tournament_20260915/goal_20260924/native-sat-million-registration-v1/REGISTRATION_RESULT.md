# Fresh SAT natural-panel registration, before dispatch

The scientific capsule was frozen from clean, pushed source commit
`9585e0370589039202901720aa2bdca1950ac35c`. Its external registration
seal is
`3b768aa18127c2aaafb261770cfafe8eed1c7631e9aaf2a42db211b678c3e470`.
The immutable source manifest SHA-256 is
`ba7c823ac3034c5d46b5c0c8e61d710da07b9e46894849c7f80b92f0e64fa1eb`.
The frozen worker SHA-256 is
`a904daf5efffa0fdb1c1e1da97e066e5729e68aabfcf29be0b46e9a7b009677d`;
the frozen controller/auditor SHA-256 is
`0719685743ffe99b0de308503fafc3451cad97b38c1dfb912599f0960f737bfd`.
The capsule is scientific (`validation_only=false`) and remains at
`/Volumes/SSD990/llm/tmp/ic-native-sat-million-registration-20261005-capsule-v1`.
Freeze built binaries and ran only `build-identity`; no ordinary query or
SAT solver was dispatched.

The complete [publication](publication-v1/PUBLICATION.json) retains the
original registration, seal, canonical config, host context and a full
6,973-file USTAR archive. The archive is 50,457,311 bytes and has SHA-256
`ecf1cd58b2e8161126b99ba1d9f6ee4541fc413212b17d2eab3820abe3a90343`.
The frozen config SHA-256 is
`42a56ed70d229b668c6b9465b5e4a54f1f82234a5e81578c846942d32bc3f70f`;
the host-context SHA-256 is
`7bd30239b7db875b552783b559877583284eb1a5229ab0ce93f0200e72b81fce`.

Both the immediate [publisher receipt](publication-local-replay.json) and
separate [data-only replay](publication-portable-replay.json) are byte-identical,
SHA-256 `d7e6940e27ab6416b8522c865393f190baff07468ba7f6e92b06863e068d6293`.
They report `PASS_DATA_ONLY_ORDINARY_SCIENTIFIC_REGISTRATION_CUSTODY`, verify
all archived bytes, run zero archived binaries and zero scientific workers,
and leave `execution_admitted=false`, `online_wall_ns=null` and
`online_speedup=null`.

This is a **before-execution** record. The protocol requires this publication,
seal and both receipts to be pushed before the one-use original frozen
controller consumes a claim. A later execution must retain every one of the
512 fixed natural attempts, including inconclusive results, and pass its own
original frozen audit. This registration makes no natural-yield, full-rank,
one-target or speedup claim.
