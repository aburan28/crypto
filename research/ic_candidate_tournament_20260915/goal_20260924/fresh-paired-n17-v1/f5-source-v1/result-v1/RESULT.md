# MatrixF5 target-free source publication

The new F5 validation capsule froze committed source
`bed4864cca41a214e8907264fb3f1d493036b80f` and the exact prior
[512-query natural F5 preparation](../../../native-ordinary-registration-v1/result-v1/RESULT.md).
That preparation has 129 verified witnesses, rank 29/29, and all 29 logs
independently replayed. The capsule's external validation seal is
`7a59a4a195426f1ece5dbdfa76849547d5258a2bae05d605d19d95a44b92f7da`.
It is validation-only, unconsumed, and has executed no target worker.

The distinct [source-only publication](publication/SOURCE_PUBLICATION.json)
has external source seal
`bd4c63599b107f97e13d26735bde87296f7800b04ac5df4120c981421df320ee`,
source-manifest SHA-256
`a34f2e9e61c78ab5ae19c75a3b0649dac6f7ae9dd7adf1c497c05da4eb5bae9f`,
worker SHA-256
`f8d64b9b9daf9097dce73d020bb196c3e73be85a765e07985913b6d33c14b749`,
and controller SHA-256
`63b3c0288d845de81e2bf414ae38ca5507eb9db05b29efb292ea1bc3edb30725`.
Its full archive SHA-256 is
`cdcd3372d4dc627aaadc630fdcf7284610f1810f3ec30ff9ff83bc8f63ff2446`.
The [independent replay](repo-replay.json) from this repository copy checked
all **6,968** archive members as data and executed zero archived binaries.
The archive contains only the immutable source, vendor, receipts, and
executables. The disclosed validation point config is absent; a small
registration sidecar retains its build identity and config digest, but no
point or seed. Appending one byte to a separate archive copy made replay
reject it as `F5 source archive bytes differ`.

The [F5 source descriptor](source-descriptor.json), SHA-256
`450f66c68d9c4df3242af1c05a35bcc1ef90f11a7470e55655e29fcb4915105d`,
is now available for a later four-arm point card. It is an assertion backed by
the complete replay above, not a fresh-target result. The frozen controller
contains a one-use card adoption path and original post-run source/card audit,
but those paths have **not** been exercised on a card. MatrixF5 still needs a
fresh card, scientific registration publication, sole target run, and
independent original audit before it can join a paired comparison. Incumbent
IC and rho also lack pre-card publications, so no card was drawn and no speedup
is claimed here.
