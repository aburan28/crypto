# Further 2× matrix-F5 work: local screen and eligible-host result

## Frozen Apple ARM64 screen

The six fixed local screens in [RECEIPTS.json](RECEIPTS.json) completed
with one release binary per screen, two seeds, one thread, three paired
repetitions per seed, and all seven F5 cases per process. All 112 processes
exited successfully and emitted seven cases. Across every arm and case,
canonical row-space fingerprint, rank, row and column counts, pruning, and
criterion work matched. Direct scalar unpack also matched raw row fingerprints,
output term counts, and reduction word operations on every call. Each
compressed receipt retains the full process output, phase times, options,
source and binary hashes, and host details. The complete exploratory source
is in [exploratory-implementation.patch.gz](exploratory-implementation.patch.gz),
generated against commit `aa38ca38eff9fcce5b51fcdb70910aea8343616e`;
the SHA-256 and byte counts are in the manifest. The submitted source keeps
only the direct-unpack option; rejected options are reproducible from the
patch without burdening the default kernel.

Numbers below are medians of three **paired complete-call** reference/new
ratios on `f5_n24_m24_d4` on an Apple ARM64 Mac. They are a nonpromoting
screen: there is no pinned x86 CPU or five-pair interval here.

| Change from current fast path | Frozen reference/new | Holdout A reference/new | Decision |
| --- | ---: | ---: | --- |
| 4 → 6 Gray-code tables | 0.971× | 0.955× | Reject; regression |
| 4 → 8 Gray-code tables | 1.001× | 0.983× | Reject; no useful gain |
| 1 → 4 pivot candidates | 0.980× | 0.985× | Reject; regression |
| 8 → 7 table bits | 1.012× | 1.010× | Reject; below 3% gate |
| 8 → 6 table bits | 1.023× | 1.018× | Reject; below 3% gate |
| Sparse leading band, dense suffix | 0.0092× | 0.0088× | Reject; 4.25 billion sparse term visits on frozen case |
| Stable lightest-first packed rows | 0.987× | 0.981× | Reject; regression |
| Direct-write scalar unpack | **1.061×** | **1.075×** | Advance to eligible-host paired test |

The sparse-leading mode reduced the frozen output from 13.73 million to
5.10 million terms, but its 8.68-second complete call was about 110 times
slower than the dense reference. This is an instructive fill result, not a
performance improvement. Direct unpack reduced the frozen local unpack phase
from about 22.3 to 18.0 ms while keeping the same 13.73 million terms.

## Eligible x86-64 one-thread run

Pending. The predeclared [paired runner](paired_f5.py) pins one allowed CPU,
requires AVX2 and BMI2, and retains all 88 processes across four seeds,
including failures. It compares the accepted opt-in fast mode with and
without `KIC_F5_UNPACK_DIRECT=1` under identical flags. The new path will be
retained only if the frozen primary paired median and lower exact five-pair
bootstrap bound exceed 1.03 and the holdout and smaller-case guards pass.
The requested further 2× remains unproved until the stronger 2.00 gate in
the [protocol](PROTOCOL.md) passes on that eligible host.

These are matrix-F5 solver-call diagnostics. They establish neither
one-target IC online time nor a DLP speedup.
