# Independent N53 C1800 striped confirmation

The striped C1800 candidate recovered and verified every scalar in the
32-run independent panel. Across eight fresh points, the C1200/striped
target S3-call geometric mean was **1.698×**, with a target-bootstrap 95%
interval of **[0.889, 3.136]**. The prespecified call-count confirmation gate
of at least 2× with a lower interval endpoint above 1 **failed**. The local
online geometric mean was **1.283×**, with interval **[0.580, 2.867]**, and
also missed the frozen exploratory 2× gate. This schedule stops after its
failed independent confirmation.

| Public point | C1200 target S3 calls | C1800 striped target S3 calls | C1200/striped calls | Exploratory online ratio |
| --- | ---: | ---: | ---: | ---: |
| T141 | 708 | 2368 | 0.299 | 0.273 |
| T142 | 1078 | 638 | 1.690 | 3.591 |
| T143 | 3710 | 1446 | 2.566 | 1.960 |
| T144 | 222 | 56 | 3.964 | 1.619 |
| T145 | 1820 | 234 | 7.778 | 10.289 |
| T146 | 890 | 886 | 1.005 | 0.974 |
| T147 | 1564 | 910 | 1.719 | 0.239 |
| T148 | 46 | 46 | 1.000 | 0.984 |

[`CONFIRM_FREEZE.json`](CONFIRM_FREEZE.json) fixes eight previously unused
public targets T141–T148, each paired on the same N53 curve and point between
the native C1200 state order and the C1800 root table with a 256-state lead
and eight striped C1200 state ranges. The candidate manifests, source and
release binary hashes are unchanged from the successful T133–T140 screen.
This protocol was committed before any confirmation timing.

Each target uses ABBA or BAAB order, two repetitions per arm, 14 preparation
workers, one target worker, a 180-second per-run cap, and a 32 GiB memory
envelope. The independent call-count confirmation gate requires all 32
scalar recoveries to verify, a target-level geometric mean C1200/striped
S3-call ratio of at least 2, and the lower endpoint of a 20,000-resample
target-bootstrap 95% interval strictly above 1. A separate exploratory
online gate requires an online ratio of at least 2 and a bootstrap lower
endpoint above 1. The local wall result remains exploratory without a
host-wide isolation receipt and passing noise gates. Preserve failures and
all target-dependent phase timings in the append-only ledger.

Run the committed protocol with:

```sh
python3 experiments/koblitz-s3-c1800-stripe-confirm-20261009/run_confirm.py \
  --binary target/release/examples/s3_partial_prepared_fullrank_c1800stripe_targetspan \
  --output /private/tmp/s3-c1800-stripe-confirm-20261009
```

The runner checks candidate and workload identities, the exact point,
full-rank preparation, factor-base and root-table identities, search
schedule, exclusive online phase sum, and independently replayed scalar.
[`preserve_confirm.py`](preserve_confirm.py) archives complete ledger rows and
raw reports. [`analyze.py`](analyze.py) checks their hashes, pairing, stage
status, and verification before regenerating
[`analysis.json`](analysis.json). The decompressed ledger SHA-256 is
`7a110869137abcae7d0e315349b25b285f26081b57fc4787c586d0340d159637`;
the raw-report archive SHA-256 is
`479605b6c4b1651ebb6efc021d5438659753cc9c6a21dedea3c64d44d60499f5`.
The local A/A 95th percentiles of absolute log online-time ratios were
2.725 for C1200 and 2.897 for striped C1800. The host record has no
exclusive partition receipt; `controlled_speedup` remains null.
