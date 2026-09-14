# Execution and recovery

The original runner recorded 198 of 288 cells, then stopped producing output for over 154 seconds. The next scheduled cell was direct S3 enumeration on fresh n18 d8 uniform target [130391, 123359]. The interruption did not localize the stall to the solver or harness and yielded no usable cost counters. An isolated replay completed six verified relations in 0.5581 seconds. This is an unresolved runner incident, not evidence of an algebraic S3 failure.

The isolated continuation runs only the remaining 90 frozen cells in the original order. It records process startup and transfer overhead separately from the three-second algorithm budget; an outer ten-second watchdog retains unknown costs as null. The affected cell remains excluded from conservative completion counts and all paired cost ratios even if the replay succeeds.

Workspace maintenance subsequently deleted the checkout, the first completed continuation and the just-started batch. Only the 198 original records had survived in a Git blob. The exact continuation and batch were rerun from published commit 833a7009299e4c4f9138934e8ca66989cf1e2827 with Python 3.12.14 and pycryptosat 5.14.7. No unavailable result is counted. Source hashes bind the recovered runs to unchanged solver files.

This is a mixed-runner cold screen with one repetition and no reserved CPU. Source/evidence preparation overlapped parts of the original cold screen. There is no timing speedup claim. The complete eight-target batch is run sequentially without concurrent timed experiments or compression/upload work; complete field primitive vectors are the main diagnostic. They still lack calibration against binary work and the full ECDLP pipeline.
