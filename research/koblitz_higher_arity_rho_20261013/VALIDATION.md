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

## End-to-end timing follow-up

| Gate | Command/surface | Outcome |
|---|---|---|
| Frozen plan | `ecbench plan --spec E2E_SPEC.json` | spec `ECS1hc7fc552628ee`; 24 executions on registered workload `W2d6c43bf8940` |
| Strict isolation attempt | `ecbench run --wait --cpus auto` | refused before execution because CPU and memory PSI were unavailable; exact stderr preserved in `records/e2e-l2-refusal.stderr.txt` |
| Diagnostic session | session `ECBS1hb77d7527b579`, explicit `--allow-busy` | 24/24 executions verified; all measured runs L0 |
| Exact audit | `ecbench verify --replay-all` | 20/20 measured deterministic replays identical; zero problems; receipt `46f230bbd9d91d8b79745a3701d547ace0f434bf2daa4232291c4d95034ffcfd` |
| Full release library | `cargo test --release --lib` | 4,059 passed, 0 failed, 114 ignored |
| Lab browser | `python3 scripts/build_lab_browser.py --check` | current after adding the new session |
| Leaderboard | `python3 scripts/build_ic_leaderboard.py --check` | current |
| Site invariants | `python3 scripts/site/test_build.py` | 48 passed |
| JSON and plot source | `jq empty ...`; `python3 -m py_compile .../plot.py` | passed |
| Benchmark copies | byte comparisons of PDF, SVG and PNG | identical to the study sources |
| Patch hygiene | `git diff --check` | clean |

The first site validation in the resumed sparse checkout could not see tracked
browser inputs. After the generator's declared scan paths were materialized,
the browser was regenerated without dropping historical identities and all 48
site tests passed. The final two-page follow-up PDF was rendered and visually
inspected page by page; labels, identifiers, tables and limits are legible.

