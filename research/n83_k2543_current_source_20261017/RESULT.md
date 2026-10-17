# K2543/root11 improves the current-source N83 exposed-point control

The current finalist IC source completed a K2543/root11 dynamic-index run on
the previously exposed N83 public point. The factor base contains **422,138**
distinct usable points and **2,543** effective columns. It reached full
rank, built **536,737,856** exact root entries, and recovered the point's
scalar from a four-point relation. Independent replay reconstructed every
base point and checked the witness sum, scalar, online phases, and exact
candidate manifest. The candidate is
`IC1N83Ckb1a1fb422138PDP4rootRCguidedLAgaussTDdirectISO0h7a345a837c55`.

| Current-source exposed-point control | K2543/root11 dynamic | K3072/root11 dynamic |
| --- | ---: | ---: |
| Actual usable base points | 422,138 | 509,952 |
| Effective columns | 2,543 | 3,072 |
| Exact root entries | 536,737,856 | 783,270,873 |
| Packed root allocation | 5,368,709,120 bytes | 5,368,709,120 bytes |
| Local cold process | 401.115 s | 678.111 s |
| One-target online | 18.983 ms | 31.908 ms |
| Independent arithmetic replay | Pass | Pass |

Both controls used source SHA-256
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`
and release binary SHA-256
`d71a9662c9d4e2dde58443bfc646a05737c6b9ca6b3d5e831b632ce09384c906`.
K2543's base hash is
`eb0a315471ae0d01598fc24d10b8dae41873cc3fab556fc91f85a584880ccf58`.
Its local online interval is 40.5% shorter than the K3072 interval on this
point; both dynamic layouts reserve the same packed root allocation. These are
shared-host exposed-point diagnostics. The next selection step compares
K2543 with K3072 and the finalist selection on previously unused public
points under alternating Tier 1 same-point IC–rho pairs.

The K2543 process had a 1,000-second cap and completed in 401.115 seconds.
It used `KIC_ROOT_ONLY11=1`, dynamic load-75 exact slots, four rank threads,
and one target thread. The producer receipt and five exclusive online phases
are versioned beside this note. The raw output, including the 33 MB base
dump, build identity, exact command and local replay scripts, is archived at
`/Volumes/SSD990/crypto/artifacts/n83_k2543_current_source_20261017/raw.tar.gz`,
SHA-256 `f0d930aed8462a951d472a0cce06a73c85e46e6fc2d1c660e073a9f09c3388e6`.
With that archive extracted to its sibling artifact directory, the
[versioned verifier](verify_archive.py) checks the archive digest and
independently replays its source, flags, entire base, relation, scalar,
phase accounting, and candidate ID:

```sh
python3 research/n83_k2543_current_source_20261017/verify_archive.py /Volumes/SSD990/crypto/artifacts/n83_k2543_current_source_20261017/raw.tar.gz /Volumes/SSD990/crypto/artifacts/n83_k2543_current_source_20261017
```
