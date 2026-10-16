# Verified official isogeny release v0.1.2390

Release: https://github.com/aburan28/crypto/releases/tag/v0.1.2390
Source: bd7eb117f094a336d3e3e8f31facb0c6b86f95bc.
Release workflow: https://github.com/aburan28/crypto/actions/runs/38915842119.
All13 Release jobs succeeded. The main CodeQL run also completed successfully.

The immutable lightweight release tag points to the exact merged source.
ReleaseMe created the draft; the byte-safe upload helper transferred the archives,
downloaded and checked every staged asset, then published. This local check
downloaded the four dedicated-tool archives afresh, matched SHA256SUMS, API
digests and sizes, and inspected safe archive members before extracting/executing.
The SHA256SUMS file itself is retained beside these receipts.

| Platform | Archive SHA256 |
|---|---|
| Linux glibc | 0dba5245855dcdedab8bf51216f5afd4ec0512847efac9a947ff6843ea0cf0a1 |
| Linux musl | 9b02dc179b360ceaeb3c6b8c349ea36fa33ad7a890ed4a67cb8a6154e0632767 |
| Apple silicon | cc88bac7575f51ea56dafeed98e40520cab5d5daac818d30e1bb4c293b6ec779 |
| Windows MSVC | 342ab3531e93e3f5cd33ec787d5438501c9f94bc0aebe3fabfd7aaa806c77ebb |

The downloaded Apple-silicon executable has SHA256
1514b1a7d7248f716e810169218dfdbdda941d454934c93e81714d68b8e66d9e.
Its version JSON equals the archive's BUILD.json, records sourcebd7eb117 and
dirty=false, and binds standard registry22b2ffef99720f2b2c7396a267412d2241a868e038a835c4fb1fd12209b102bf.
The embedded challenge policy exposes54 standard-prime boards. The CLI package
version remains0.1.0; the rolling repository release tag isv0.1.2390, with exact
source revision providing the executable identity.

From a fresh directory outside the source checkout, the downloaded native
binary admitted the existing P192/13 and P521/7 public control artifacts.
All20 scalar transports per map passed, and both supervisor receipts are PASS
with exit0 and zero unavailable RSS samples. No new isogeny search, large
target, benchmark or scientific result was run. Linux/Windows archives were
locally checked for integrity; their runtime checks are the successful CI jobs,
not execution on this Apple host.

For Apple silicon, extract isogeny-v0.1.2390-aarch64-apple-darwin.tar.gz and run:

    ./isogeny-v0.1.2390-aarch64-apple-darwin/isogeny --version
    ./isogeny-v0.1.2390-aarch64-apple-darwin/isogeny curves
    ./isogeny-v0.1.2390-aarch64-apple-darwin/isogeny --help

No Cargo, checkout, external registry or companion executable is required.
The default-off fixed-modulus optimization is still in PR1991; it is not part
of release2390. Its private matched binaries and S3 custody question remain
separate from this already authorized public release.

The separate reports-s3 workflow failed in credential configuration before
report upload with an AWS signature mismatch. Its exact log is retained;
reports-s3.yml is unchanged by the repair. This is not a green report-mirroring
outcome, and credentials were not changed.

raw.tar.gz retains40 ordinary evidence files,323238 bytes, SHA256
d9eed55fde97df3d05b0ab0946f268a2690368d2f95003af41f385d3cc5a33aa.
All listed hashes replayed. Executable/archive payloads remain outside Git;
the immutable GitHub asset IDs and hashes are in download-integrity.json.
The source-linked diagram and inspected report supplement these raw receipts.
The old v0.1.2372 bytes and warning remain preserved.


Exact public source identities:

| Control | ICV1 source |
|---|---|
| p192 | icv1-fp192-t31607402316713927207482677199-52e4af59 |
| p521 | icv1-fp521-t657877501894328237357444332315020117536923257219387276263472201219398408051703-f767e0d1 |
