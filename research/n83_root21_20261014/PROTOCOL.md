# N83 exact 21-byte root-pair capacity gate

Two N83 normal-basis roots contain 166 bits. The new optional
`KIC_ROOT_ONLY21=1` layout concatenates them in 21 bytes and recovers the
state labels from the deterministic symmetric traversal ordinal. The indexed
five-byte root table can request exactly `2^29` entries, yielding `2^30`
slots, when `KIC_INDEXED5_FIXED30=1`. Construction checks that twice the
complete state bound is below the slot count. Full-key equality still reads
the retained roots, and first-inserted witnesses remain the lookup policy.
For indexed insertion, the source now obtains each root directly; ordinal
label recovery occurs only when a matched witness needs it.

## Frozen instance and variants

- Curve: `y^2+xy=x^3+x^2+1` over `F_2^83`, polynomial low terms
  `[0,2,4,7]`; subgroup order `8569786107849059`, cofactor `1128547018`,
  generator `(6614244651398184860292015,8868260410029249391539529)`.
- Public point: `research/n83_indexed5_20261013/primary_target.jsonl`,
  SHA-256 `21c8dfc7d488909e5c4a42c839fcc01bdefa9821c11ba001626a1d15697240e8`.
  The scalar sidecar SHA-256 is
  `651eb7d0ef08a1655036746311b19b6989e340ca4752b6d675606055e59c12a1`.
- Producer source SHA-256:
  `accd99313204d0e3327968230ed3c64df49d661ebddb466147a9a8fc844bd1c9`.
  `Cargo.lock` SHA-256:
  `ac5cd7a3459be776180afe1eb2975c76b1a88f84bd1e99aa20475927ec260a90`.
- All arms use symmetric S3, exact five-byte root slots, no Bloom filter,
  no sparse canonicalization, rank seed 29, four rank threads, one target
  thread, and a fresh process. The first four arms use 21-byte roots and
  fixed `2^30` slots; the last uses 22-byte roots and the existing
  `ceil(4S/3)` request. The two K2543 arms are exact-source controls.

| Order | Arm | K | Predicted usable points | Complete state bound | Array payload GiB | Child cap |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 1 | root21 fixed30 | 3240 | 537,840 | 435,652,020 | 13.520 | 2400 s |
| 2 | root21 fixed30 | 3220 | 534,520 | 430,290,210 | 13.416 | 2400 s |
| 3 | root21 fixed30 | 3190 | 529,540 | 422,309,745 | 13.259 | 2100 s |
| 4 | root21 fixed30 | 2543 | 422,138 | 268,375,505 | 10.248 | 900 s |
| 5 | root22 load75 | 2543 | 422,138 | 268,375,505 | 10.499 | 900 s |

The bound is `83*K*(K-1)/2+42*K`. Lean 4.34.1 proves exact root-pair
roundtrip, that both roots fit in 21 bytes, and that K3240's complete key
upper bound and state IDs fit the fixed table. It proves the K3240 array
payload is below 13.6 GiB; allocator and other process memory are measured
separately. The completed K3114 22-byte cloud arm measured a 13.637-GiB
`wait4` peak against its 13.245-GiB array payload. The runner stops a child
above 14 GiB sampled RSS and records the final `wait4` peak even if it
exceeds the cap between samples. Every arm, including failures and later
fallbacks, retains stdout, stderr, usage, cgroup, and wall-cap receipts.

The same-point strong-rho arm uses archived `fold-L64-DP6-E0` config SHA-256
`cf76ebb5372acfc1355b88374f875980be66cdffba1a82cb6e4c2f28bcfa723b`,
seed `2026101301`, one thread, and a 120-second cap. The independent
verifier checks each completed base, subgroup, four-point relation, scalar,
exact state count, allocation and five exclusive online phases. It checks
K2543 root-key count, base hash and scalar against the 22-byte control.

The one-run cloud-instance timing is exploratory. After the largest
in-envelope base is selected using the archived secondary target panel, a
new public point and five alternating pairs will test the selected
candidate against same-point rho under the Tier 1 timing gate.
