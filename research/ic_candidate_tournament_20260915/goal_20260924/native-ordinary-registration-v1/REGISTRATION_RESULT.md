# Scientific ordinary-panel registrations: data custody, before dispatch

Both arms have been frozen from source commit
`4e8634896cfa6271c550013ab99415e5bc771a02` and the common source
manifest SHA-256
`61a1f8e33b30e4a1234490acfde4fa482622779d3037e3bf639963a8f3735375`.
The scientific freeze (`validation_only=false`) built a release worker and
original checker, retained offline dependencies and build receipts, then
verified its immutable inventory. It ran `build-identity`, not an IC query.
The original local capsules remain at
`/Volumes/SSD990/llm/tmp/ic-native-ordinary-registration-20261005-{f5,cms}-capsule-v1`;
the complete portable archives and sidecars below are the durable PR evidence.

| Arm | Frozen config SHA-256 | External registration seal | Full archive SHA-256 | Bytes | Regular files |
| --- | --- | --- | --- | ---: | ---: |
| MatrixF5 | `84c13f310c03aa72616f5926453950bfabdb31a8e00a229bf331b2a67d435e6f` | `89de6a771df622536b83ae15744f47a429c3e63063a696c6b9012e8ca0fcca97` | `010b4ceb0db720168f90b99166ce2be9b5e5a7a5c4233f94210b2f9faa235293` | 35,625,082 | 6,948 |
| External CryptoMiniSat | `9313ffd80badfb432b4bc466786255d6afc69555b5a589ed8d3745df861bfead` | `8da6697ce17967648a4f72ef9d8bd4020cf3f877575369dc0e5afdfad95b20bf` | `ac9f4c284f14530c8c5e3629607b05cf9ed7c28ba2bad279f364211f7e841bfc` | 50,224,359 | 6,965 |

The actual [F5 publication](f5-publication-v1/PUBLICATION.json) and
[CryptoMiniSat publication](cms-publication-v1/PUBLICATION.json) retain their
`capsule.tar.gz`, config, host context, registration, seal and build receipts.
Each archive contains every sealed immutable file, including the exact source,
vendored dependencies, compiled binaries, and the accepted CMS archive where
applicable. The [F5 publisher receipt](f5-publication-replay-v1.json) and
[portable replay receipt](f5-portable-replay-v1.json) are byte-identical; the
[CMS publisher receipt](cms-publication-replay-v1.json) and
[portable replay receipt](cms-portable-replay-v1.json) are likewise
byte-identical. Each replay verifies the full USTAR inventory as data without
extracting or executing an archived binary.

The host-context SHA-256 is
`7bd30239b7db875b552783b559877583284eb1a5229ab0ce93f0200e72b81fce`.
Both receipts say `PASS_DATA_ONLY_ORDINARY_SCIENTIFIC_REGISTRATION_CUSTODY`,
`scientific_worker_calls=0`, `execution_admitted=false`, and
`online_wall_ns=null`. No natural-yield, full-rank, target-solve, rho or
speedup claim follows from registration. The original capsules have not been
consumed. The validation-only archive is a separate historical control and
scientific replay rejects it.

The two one-use executions may start only after this publication, seals,
configs and source commit have been committed and pushed. The original frozen
checker must receive `--publication` and its exact external seal. Keep its
terminal files and an audit made by that same original checker. The
[protocol](PROTOCOL.md) fixes the stop, outcome, accounting and later
complete-solve gates; a failed or partial arm is retained without retry.
