# Validation log

## Final passing gates

| Gate | Command/surface | Outcome |
|---|---|---|
| Focused release library tests | `cargo test --release --lib koblitz_sixsum_extract -- --nocapture` | 6 passed, 0 failed, 4,167 filtered |
| Ecbench plan | `ecbench plan --spec SPEC.json` | 20 executions; all five arms resolved on the registered n=41 workload |
| Ecbench session | session `ECBS1h4f09e1612522` | 20/20 verified |
| Exact audit | `ecbench verify --replay-all --exit-code` | 15/15 measured deterministic replays identical; zero problems |
| Leaderboard | `scripts/build_ic_leaderboard.py --check` | current |
| Lab browser | `scripts/build_lab_browser.py --check` | current |
| Site invariants | `scripts/site/test_build.py` | 48 passed |
| Patch hygiene | `git diff --check` | clean |

The focused release tests were run on the final candidate logic. A later
format-only reduction restored unrelated pre-existing line wrapping without
changing the Rust AST or measured path.

## Preserved failures and downgrades

1. `conductor check` could not run because this runtime has no `conductor`
   executable (exit 127); no `.conductor/WORKFLOW.md` exists.
2. The first release-test compilation failed because the partial clone omitted
   pre-existing `include_str!` research fixtures. Only the exact named fixture
   directories were materialized; the same test then compiled and passed.
3. A proposed infinity-branch direct root was tested after the first successful
   session. The focused test rejected it: candidate roots
   `[72577, 105885, 126501]` versus affine quartic roots
   `[72577, 126501]`. The proposal was reverted, and the final 6/6 release
   suite passed. This confirms that the original infinity skip matches the
   affine S4 polynomial.
4. Ecbench first refused the session because PSI CPU and memory telemetry are
   unavailable. The explicit `--allow-busy` rerun earned L0; counts remain
   primary and wall results are descriptive.
5. The n=83 factor-base gate refused with
   `no factor-base plug-in runs on a wide (n > 62) curve yet`. No n=83 IC
   timing or scaling value was invented.
6. SVG rendering through ImageMagick and Chromium was unavailable in this
   sandbox. The committed PNG was rendered from the equivalent TeX chart
   source, and the two-page report PDF was visually inspected after rendering.

