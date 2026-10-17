# Validation scope

- Native archival auditor:16 tests passed;112 payloads and190 recorded context-output comparisons verified. No Python research execution or new timings.
- Four PDFs / seven pages rendered and visually inspected; benchmark copies byte-identical to study originals.
- Root `cargo test --release --lib --offline`: 4,125 passed, 10 failed, 117 ignored. One failure was a sparse-checkout fixture; nine were local-socket permissions in the restricted sandbox. After restoring tracked fixtures and enabling local sockets, all 37 pollard_collab tests and the one affected fixture test passed. No test/source was changed or skipped to achieve this. The original failed full run is retained; no single clean full-suite pass is claimed. No Rust library or example implementation was touched.
- Three existing scoreboard-only test methods passed; generated panel index is current at303panels. Whole-site assembly was not run.
- A first scoreboard test-reporter invocation used the wrong stream wrapper; the preparation error is retained separately, followed by the three passing checks.
- Historical oracle records are preserved; the former Python suites were not rerun as new research, consistent with the native-only rule.
- Conductor CLI and tracked .conductor/WORKFLOW.md unavailable at the pinned base. No task reservation is claimed.
- No independent equation re-solve, qualified wall result, matched F4/F5/F6 comparison or scientific promotion.
