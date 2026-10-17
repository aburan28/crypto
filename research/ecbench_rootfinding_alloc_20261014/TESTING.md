# Verification and integration receipts

- `cargo test --offline --locked --release --lib prime_lane_collector_has_identical_serial_and_batched_streams -- --test-threads=1` passed before the candidate binary build and again after merge commit `16d0cd44f`; the postmerge raw output is [final-lane-test.log](final-lane-test.log.gz). It checks ordered points and coefficients over 2,048 probes at each of widths 4, 8, and 16, including restarts.
- The collector source blob was `7d078bcfad692437db8ab159c7a637861a7c7b99` both at evidence commit `0ec9e598f` and after merge commit `16d0cd44f`, so the final focused test covers the same implementation as the sealed candidate binary.
- The baseline and candidate plans compare byte-for-byte. Each of the three sessions has 128/128 verified records and a successful `verify --replay-all --exit-code` receipt with 96/96 measured replays identical; see [REPORT.md](REPORT.md).
- `python3 compare_sessions.py` on candidate/baseline and repeated-baseline/baseline returns zero mismatches in 128 compared deterministic records per comparison. The initial checker error and corrected rerun remain separate.
- `python3 scripts/build_lab_browser.py --check` passed and printed `docs/browser/data.json is current`; [log](browser-check.log).
- `python3 scripts/build_ic_leaderboard.py --check` exited 0 without output; [empty log](leaderboard-check.log).
- `python3 scripts/site/test_build.py` passed 48 tests; [log](site-test.log).
- `typst compile` produced the combined [PDF](REPORT.pdf) and the three session PDFs. The combined PDF was visually inspected after rendering to PNG.

The failed first release build is recorded in [BUILD_FAILURE.md](BUILD_FAILURE.md); the retry with `TMPDIR` on `/Volumes/SSD990-2` succeeded. The old baseline binary and new candidate binary SHA-256 values are bound in the [report](REPORT.md) and session files. All timed runs preceded the postmerge site and source checks.
