# Stage 27 direct-MITM results

Stage 27 completed the matched direct meet-in-the-middle lane for all 160 truth-free Stage 26 inputs. GitHub Actions run [34633920325](https://github.com/aburan28/crypto/actions/runs/34633920325) executed from merged commit `02225ed61e7331c11a1e4665b953a8384b0b08de`. Its four cells ran concurrently, while every cell restricted its packet verification and decomposition process tree to one Linux CPU.

The post-seal scorer authenticated every copied source against packet inventory `c937afb9b172d114768b0b96a4b5ccf66e91fbf278b77c58ed49f0fd74af37f7`, reconstructed all 160 task records, checked every exhaustive terminal and exact three-point SAT witness, and only then opened the committed Phase A oracle. Each cell produced 20 true positives and 20 true negatives. There were zero false or inconclusive outcomes.

| Cell | Factor points | Direct process core-s | Inclusive one-CPU elapsed-s | Peak process RSS | Peak sampled tree RSS | Group additions |
|---|---:|---:|---:|---:|---:|---:|
| n31 standard | 31 | 0.643425 | 187.936816 | 14,393,344 | 51,269,632 | 20,645 |
| n31 GGMP | 63 | 0.703140 | 140.530319 | 14,393,344 | 51,900,416 | 82,178 |
| n41 standard | 25 | 0.629727 | 144.022971 | 14,393,344 | 52,867,072 | 13,632 |
| n59 standard | 483 | 161.336348 | 294.048550 | 38,166,528 | 85,479,424 | 4,687,797 |

The direct-MITM children consumed 163.312640 core-seconds. Inclusive cell envelopes consumed 759.710741 core-seconds and 766.538655 summed one-CPU elapsed-seconds because each independently charged packet verification, setup, source verification, factor-point materialisation, pair-table construction, and target search. Parallel cell job span was 307 seconds; complete workflow wall time was 322 seconds. The panel constructed 4,768,800 pair-table entries and charged 4,804,252 group additions. SAT conflict counts do not apply to this algorithm.

The backend binary came from the immutable tool artifact built in Stage 26 run 34632018379 and was authenticated by SHA-256 `667fc00681e114acaaf662a8c00c214354e1a08e236a6974fb08e52504711ea0`; its build cost is therefore charged once rather than duplicated. The four downloaded artifacts are retained in the terminal-evidence archive with SHA-256 `5a48ef9004f896e3a9eee16a170bbd9b135accd9fc6aa9fe002f38602ffd777a`.

This closes the same-instance direct-decomposition lane. It remains a finite public PDP result, not a complete end-to-end index-calculus cost or a Koblitz SOTA result. The complete licensed Magma panel, unified full-cost composition, and unaffiliated reproduction and novelty review remain open.
