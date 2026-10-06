# Post-buffer-marked CryptoMiniSat build

The guarded second attempt built a new `prepared-cms` executable after the
original one-use SAT preparation audit passed at 29/29. This is a build
receipt, **not** a transport-parity pass, target solve or speed measurement.
The first attempt and its corrected archive-hash transcription are recorded
in [BUILD_ATTEMPTS.md](../BUILD_ATTEMPTS.md).

The binary SHA-256 is
`128f074850bbd9d973bd50be26c6ffc161e97c7a365b9234f8af171545fe7263`.
It is distinct from the accepted cold CMS binary. The original build tree is
retained at `/Volumes/SSD990/llm/tmp/ic-marked-cms-build-20261006-v3`; the
portable executable and small original build logs/receipts are copied here
byte for byte. `build-inputs.json` pins the accepted source, CaDiCaL and
CaDiBack archives, the post-buffer patch, compiler binaries, preparation
audit, terminal and corrected build-script SHA-256. The original accepted
source/dependency bytes remain in the archived native-assets bundle under
`native-sat-million-registration-v1/publication-v1`; the exact extraction,
patch and compile commands are in [build-marked-cms.sh](../build-marked-cms.sh).

The terminal status is `BUILT_UNVALIDATED`, with transport parity and target
execution admission false. Before this binary is eligible for a target
capsule, run the predeclared nine disclosed roles with their fixed 120-second
watchdogs, then independently audit the original files as data. Preserve a
failure or timeout as a failed gate; do not substitute a diagnostic rerun.
