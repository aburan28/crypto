# Sparse canonical selection and Bloom compose on N83

With the exact K3300/root11 factor base fixed, adding the three-hash Bloom
prefilter to sparse canonical selection reduced exploratory online time by
**3.644 ms per target** on 1,000 paired public subgroup points. The frozen
contiguous-block bootstrap 95% interval for Bloom-on minus Bloom-off is
**[−3.834, −3.447] ms**. Both arms independently replayed all target answers
and four-point relations, and both remained below the 14 GiB RSS guard.

| Configuration | Mean online time per point | Preparation | Peak RSS |
| --- | ---: | ---: | ---: |
| Sparse canonical, Bloom off | 32.572 ms | 467.714 s | 10.05 GiB |
| Sparse canonical, Bloom on | 28.928 ms | 487.654 s | 10.55 GiB |

The pair uses the same 547,800 factor-base points, 3,300 folded columns,
903,852,075 exact root keys, source snapshot, and ordered target list. Bloom
added 19.940 seconds of target-independent preparation and 0.50 GiB of peak
RSS in this run. The online gain is 11.19% relative to sparse without Bloom.
Each arm started from an empty index and log cache. The child execution cap
was 1,200 seconds; the observed peak RSS values were 10,792,984,576 and
11,329,900,544 bytes.

The selected combined candidate is
`IC1N83Ckb1a1fb547800PDP4rootRCguidedLAgaussTDdirectISO0h2aa3ee55ba07`.
This panel is a secondary, reused-point cloud-instance diagnostic. The
five fresh, alternating same-point Tier 1 IC/rho pairs are recorded in the
separate Bloom reserve study when its run completes.

The producer was GitHub Actions run `38679320666`, artifact `11948835612`.
Its archived ZIP SHA-256 is
`20be795e5e1a596f92f357e8764d2529ee349a546d18c439b6bd8354901d1916`;
the independent analysis SHA-256 is
`cb5ba82cab460faa12b27c5fe8a7926899be538d18fb4d81c44bb8760f71267f`.
The executed source commit is `9289963251767f88340e1c8900d18f50f1b40ee0`,
the IC source SHA-256 is
`814afa03ae1c386e62e518be3c7a55fb108b09c41b05bef64ca1371787466c39`,
and the shared workload ID is `0d7b05d9fca6`. The raw ZIP is retained at
`/Volumes/SSD990-2/n83-archive-20261015/raw/n83-root11-bloom-sparse-38679320666.zip`.
The candidate manifests, compressed per-target measurement rows, selection,
and record index beside this report bind the result to those receipts.
