# Fresh-seed Linux confirmation: PASS

The frozen [protocol](PROTOCOL.md) passed on [Actions run 37148557781](https://github.com/aburan28/crypto/actions/runs/37148557781) at PR #1292 head `7151aea471f824f60fd4b9c60d3341c4ccfdd46c` (Actions merge SHA `b0421000d40e7657afbbb0cbc2f322056ffa0be8`). The runner was native x86-64 Linux on AMD EPYC 7763 with AVX2/BMI2, and rebuilt with Rust 1.98.1. Release GF(2) and F5 rank tests passed. These are Boolean matrix-F5 solver-stage **complete-call** timings; they do not measure IC online time or a rho comparison.

The reference was selective echelon; the candidate was the exact cut-19 support certificate with width-8 row-basis tables. Each cell used one warmup per arm, five A/A pairs and five alternating paired calls. Reference and candidate wall columns below are medians of their five paired calls in milliseconds; speedup is the median of five *paired ratios*. The lower bound is the exact five-pair bootstrap 95% lower endpoint. All figures exclude process launch and match the frozen call boundary.

| Threads | Fresh seed | Selected attempt | Reference ms | Candidate ms | Paired speedup | 95% lower | Candidate XOR / reference |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | `confirm_a` | 1 | 104.521 | 34.491 | 3.019× | 2.999× | 28.095% |
| 1 | `confirm_b` | 2 | 105.054 | 34.269 | 3.059× | 2.995× | 28.095% |
| 1 | `confirm_c` | 2 | 105.777 | 34.458 | 3.070× | 3.019× | 28.078% |
| 1 | `confirm_d` | 2 | 103.813 | 34.404 | 3.018× | 3.006× | 28.078% |
| 2 | `confirm_a` | 1 | 95.474 | 35.125 | 2.720× | 2.703× | 28.095% |
| 2 | `confirm_b` | 1 | 96.208 | 34.816 | 2.763× | 2.700× | 28.095% |
| 2 | `confirm_c` | 1 | 95.335 | 34.965 | 2.748× | 2.672× | 28.078% |
| 2 | `confirm_d` | 1 | 95.506 | 34.994 | 2.751× | 2.650× | 28.078% |

Every selected attempt was complete and isolated with zero contended samples and no foreign user threads on reserved CPUs. The first one-thread attempts for `confirm_b`, `confirm_c` and `confirm_d` were rejected by the isolation preflight because of high pressure stall information; they made zero benchmark calls and are retained in the artifact. No seed or attempt was selected by its timing. Both native analyzer reports say `pass` for all eight cells, with no errors. They checked the fixed cut, table width, source hashes, rank, canonical row space, criterion, built/pruned counts, and original-row certificate route. All six smaller cases on all 22 calls had the certificate route inactive and equal complete non-timing output; their A/A and paired wall ratios remain in the raw receipts as diagnostics.

The full compressed artifact, with every attempted run, CPU choice, isolation receipt, 22 per-call outputs per selected cell, source and binary hashes, and both analyzer reports, is `linux_run_37148557781.tar.gz` (SHA-256 `fd2d5a25bdef344584b2f5a71f0c857b239b1fbf771278bb0f129f1017a45478`). Extract with `tar -xzf linux_run_37148557781.tar.gz -C DIR`. The measured F5 source SHA-256 is `7a0b16799932b5f3d814621c989a5e15afcca8a2fd1feef91d1ac1dd09c9b422`; the GF(2) kernel source SHA-256 is `4b5cacbcc7f1150ed950d5d468c6847678f73fa3a93b736863ebc6f5a0dd2bf9`.

The earlier four-seed width-8 and width-6 experiments **failed their own original overall gates** on smaller-case wall controls. Their failures and raw receipts remain at `../f5_support_cut19_call_20261002/LINUX_RESULT.md` and `../f5_support_rank_bits_20261002/LINUX_RESULT.md`. This fresh, independently seeded confirmation passes the revised preregistered gate; it does not relabel those failures.
