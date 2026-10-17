#set page(paper: "us-letter", margin: 0.65in)
#set text(size: 9.6pt)
#set par(leading: 0.43em)
= Verified native isogeny release v0.1.2390
Source bd7eb117f094a336d3e3e8f31facb0c6b86f95bc. All 13 Release jobs and the main CodeQL workflow passed. The published non-draft tag points to that exact source.

#image("publication.svg", width: 100%)

Fresh downloads of all four dedicated-tool archives match SHA256SUMS, GitHub asset digests and API sizes. The local checker inspected archive members before extracting or executing. Full hashes, asset IDs and file lists are in download-integrity.json.

#table(columns: (1.3fr, 2.8fr), inset: 5pt, stroke: 0.4pt + rgb("#b7c6d0"),
 [*Gate*], [*Observed result*],
 [Exact release tag], [bd7eb117 source; draft=false; v0.1.2390],
 [Release workflow], [13 successful jobs; native platform smoke tests included],
 [Fresh archive integrity], [Linux glibc, Linux musl, Apple silicon, Windows: all 4 match],
 [Downloaded Apple version], [BUILD.json matches runtime version; dirty=false; registry22b2ffef],
 [Fresh native admission], [P192/13 and P521/7 ACCEPTED; all 20 transports per map],
 [Memory-monitor receipts], [Both PASS / exit0 / zero unavailable samples],
 [Standalone use], [Fresh cwd; no Cargo, checkout or external registry],
 [Separate report mirror], [AWS credential-signature failure before upload; retained],
)

The Apple executable SHA256 is
#text(size: 8pt)[#raw("1514b1a7d7248f716e810169218dfdbdda941d454934c93e81714d68b8e66d9e")].
Its archive SHA256 is
#text(size: 8pt)[#raw("cc88bac7575f51ea56dafeed98e40520cab5d5daac818d30e1bb4c293b6ec779")].

*Run locally.* Extract the matching dedicated archive and invoke #raw("isogeny --version"), #raw("isogeny curves") or #raw("isogeny --help"). Runtime package metadata remains 0.1.0; the rolling release tag and exact Git commit provide the release identity.

*Scope.* Release 2390 contains the merged compact-certificate, tracing and byte-safe publication work. The experimental fixed-modulus feature remains in PR1991 and is not in this release. Its four local research executables are not uploaded by this QA step. No new isogeny, benchmark or board promotion was run.

*Evidence.* raw.tar.gz contains 40 hash-replayed ordinary files;323238 bytes, SHA256
#text(size: 7.7pt)[#raw("d9eed55fde97df3d05b0ab0946f268a2690368d2f95003af41f385d3cc5a33aa")].
The large archive/executable payloads stay outside Git. Local runtime verification covers Apple silicon; Linux/Windows runtime evidence is the successful CI jobs. Old failed assets and the reports-s3 credential failure remain preserved.

#link("https://github.com/aburan28/crypto/releases/tag/v0.1.2390")[Official release] · #link("https://github.com/aburan28/crypto/actions/runs/38915842119")[Exact Release run]

