# K1101 capacity-edge producer and independent replay

The merged N83 producer completed a one-point cold process at K1101 on the
local macOS host in **107.760 s** of producer `process_total`. The stage
intervals were 42.722 s setup, 64.733 s relation/rank and 0.300 s target.
The separate outer child wall was 108.122 s. This is a correctness and
layout diagnostic on an ordinary host; it is not a Tier 1 factor-base
comparison.

The exact source SHA-256 was
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`,
the compiled macOS binary SHA-256 was
`7b91925a62af09a6fa7f3393bc12c252ad0f38093125df42883f2bd64881fef2`,
and the input was the already exposed public point frozen in the
[grid protocol](../n83_e2e_factor_base_20261016/PROTOCOL.md). The producer
enumerated **182,766** subgroup-usable points, attained rank **1101**,
allocated **671,088,640 bytes** for the exact root table, and recorded base
hash `20676b872a392ee551fb8a300896b64fbdab56565a77240ebf5a61dbbffb7dd7`.
Independent group arithmetic replayed the enumerated base, four selected
point indices, their sum, the recovered subgroup scalar and the five
exclusive target-phase costs. The replay status was
`PASS_INDEPENDENT_BASE_RELATION_SCALAR_REPLAY`.

The complete raw output, base dump, command flags, stderr, binary and
independent replay receipt are archived at
`/Volumes/SSD990-2/n83-archive-20261016/raw/n83-k1101-local-smoke-20261016.tar.zst`,
SHA-256
`40ed8bf8cdd8b554f7fc8a7d547bc47f47dde242ae4dcf6cd6e57ac00badfe27`.
The `zstd` integrity check and complete tar listing passed. The replay JSON
SHA-256 is
`edc0c6be34d1098cd5e43ec8f6c3215df6e534567911cc1f051662fcb5681275`.
After extracting the archive, rerun the independent checker:

```sh
python3 research/n83_e2e_holdout_20261016/verify_local_k1101.py \
  /path/to/n83-k1101-current-source-smoke-20261016
```

It replays the base and four-point relation and checks the archived receipt
exactly.
The 20-value Tier 1 cloud screen and five fresh-point holdout remain the
selection gates.
