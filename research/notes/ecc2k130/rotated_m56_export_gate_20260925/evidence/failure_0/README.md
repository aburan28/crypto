# Preserved first frozen gate attempt

Command: `python3.12 research/notes/ecc2k130/rotated_m56_export_gate_20260925/run.py --out /private/tmp/kic-rotated-m56-run-20260925` on preregistered head `572ebd1af351ccd09667a674f09c4f994b25a9ff`.
The runner started the n13 child at 2026-09-25 12:40:53.075810 UTC and it exited 1 at 12:40:53.133707 UTC, before either tuple panel ran. The n19 child was not started. `receipt.json`, the exact child JSON, stdout and stderr are copied byte-for-byte from that failed directory; its own per-file SHA-256 values and source command are in the receipt. The n13 stderr identifies a preflight assertion in `gate.py` line 315.

The assertion compared the factor list's x-values with **every** x-value in its two-dimensional subspace. Some x-values do not lift to rational points, so #767's correct factor list is a proper subset. The subsequent source amendment retains the independent complete-factor-list comparison, verifies subset membership, and counts only masks with rational lifts. This attempt yielded no tuple, trace, chain, symmetry, solver or performance outcome. It remains preserved and is not overwritten by the amended run.
