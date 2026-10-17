# Three adjacent N83 root11 bases pass full single-point arithmetic replay

K2752, K2880 and K2944/root11 each constructed an exact base, reached full
rank, recovered the scalar of the same previously exposed public point, and
passed independent replay of every base point, its four-point relation and
the scalar. This admits all three to the separately frozen 128-point
secondary query-distribution batch. Their exact current-source identities
and local feasibility measurements are:

| K | Usable base points | Folded columns | One-target online, ms | Complete process, s | Packed root state, GiB | Exact IC1 candidate ID |
| ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 2752 | 456,832 | 2,752 | 7.113 | 872.833 | 3.220 | `IC1N83Ckb1a1fb456832PDP4rootRCguidedLAgaussTDdirectISO0hcbd26a90fa18` |
| 2880 | 478,080 | 2,880 | 16.792 | 537.680 | 3.526 | `IC1N83Ckb1a1fb478080PDP4rootRCguidedLAgaussTDdirectISO0h8daf9a7cfc40` |
| 2944 | 488,704 | 2,944 | 21.694 | 657.156 | 3.685 | `IC1N83Ckb1a1fb488704PDP4rootRCguidedLAgaussTDdirectISO0h042446828aa7` |

The three `screen_result_k*.json` files, exact candidate manifests in the
parent directory, and `candidate_manifest_receipt.json` bind each result to
the base hash and raw command. `plan.json` fixed the sizes, the 900-second
cap, executed source SHA-256
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`
and binary SHA-256
`faaac22aa6228c8d486875d6947fa654417d1122191fbe0eb9bde59bb45eb98f`
before the runs. `screen.py` and `materialize_candidates.py` are the exact
executed source snapshots; their SHA-256 digests match the respective frozen
plans. The raw base dumps, target records, command environments and stdout
remain under
`/Volumes/SSD990/crypto/artifacts/n83_root11_finegrid_20261017/`.

These timings come from one exposed point on an ordinary shared macOS ARM64
host with 48 GiB RAM. K2752's index-build phase alone took 825.3 seconds,
versus 404.8 seconds in an earlier K2816 local run under the same algorithm
flags; the host was not isolated, so the process times are feasibility
records. The frozen fresh 128-point batch will compare online query
distributions with one shared setup per arm. A later Tier 1 holdout will
measure one-point IC/rho pairs with independent host calibration.
