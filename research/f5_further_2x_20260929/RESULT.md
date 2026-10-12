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

[CI run 36601441800](https://github.com/aburan28/crypto/actions/runs/36601441800)
completed all 88 same-binary processes on one pinned CPU of an AMD EPYC 7763
with AVX2 and BMI2. The release F5 tests passed. Every process succeeded;
on all seven cases and four seeds, prior and new matched raw and canonical
row fingerprints, rank, output term count, row and column counts, pruning,
criterion work, and reduction word operations. The complete 804,971-byte
receipt is [runs/36601441800-t1.json.gz](runs/36601441800-t1.json.gz),
with compressed SHA-256
`594cfb3bf0cd4339119b0b3df824a7a3267f5243c7e1d15bc1c5436c9e0c5e10`
and uncompressed SHA-256
`6389462398f461c677179840d0508298f1f8aa936b3031323f1ea09dadc96bca`.
It retains every process output and status, source and binary hashes, CPU
features, affinity, load, phase timings and signatures. The Actions artifact
retains the same full JSON. The measured source SHA-256
`bf35680f431ab7ac13a7615d334ec9b060217b84bdde36f03ce00319491d42da`
matches the submitted `matrix_f5_f2.rs`.

The reference and candidate both select selective echelon output, fused row
counting, direct packed rows, AVX2 XOR, and table reuse. Only the candidate
sets `KIC_F5_UNPACK_DIRECT=1`. The table reports medians of five paired
**complete-call** ratios, exact five-pair bootstrap 95% intervals, and each
seed's reference/reference A/A maximum. Larger ratios are better.

| Seed | Prior / direct unpack (95% interval) | A/A maximum |
| --- | ---: | ---: |
| Frozen | **1.048× (1.031–1.056×)** | 1.014× |
| Holdout A | 1.042× (1.034–1.075×) | 1.004× |
| Holdout B | 1.041× (1.008–1.045×) | 1.012× |
| Holdout C | 1.040× (1.034–1.048×) | 1.005× |

On the frozen primary, marginal median complete calls were 128.31 ms prior
and 122.96 ms with direct unpack on this runner. The paired unpack-phase
ratio was 1.099× (1.091–1.110×). All smaller-case complete-call medians
exceeded their own A/A minima; the smallest margin was 0.0053× on frozen
`f5_n12_m12_d4`. The predeclared **incremental** gate passes, so the path
remains opt-in. The further 2× gate fails: the frozen complete-call ratio is
1.048×, far below 2.00. Reaching that gate on this workload would require
about 64 ms or less against the 128.31 ms reference under matched resources.

The one-off CI workflow source is archived as [WORKFLOW.yml](WORKFLOW.yml)
after its successful run; it no longer triggers on every later F5 PR.

These are matrix-F5 solver-call diagnostics. They establish neither
one-target IC online time nor a DLP speedup.
