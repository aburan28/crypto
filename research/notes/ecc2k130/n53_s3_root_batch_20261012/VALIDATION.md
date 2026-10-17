# Validation and known baseline site error

All commands ran from the repository root on 2026-10-12 unless a subdirectory
is specified. The valid experiment execution used the frozen binary/source
hashes in [`FROZEN.json`](FROZEN.json), and the native analyzer generated
[`ANALYSIS.json`](ANALYSIS.json). The separate sandbox-wrapper failure is
recorded in its [run report](reports/run-20261012T1735Z/REPORT.md).

| Command or check | Result |
|:--|:--|
| `cargo test --release --lib` | 4,000 passed, 0 failed, 113 ignored; 138.42 s on final branch |
| `cargo test --release --example koblitz_s3_root_batch_experiment --example n53_s3_root_batch_run --example n53_s3_root_batch_report --example koblitz_rho_fixture --example koblitz_rho_batch_ks_strong_online --example koblitz_one_target_replay` | 22 passed across six release example targets, 0 failed |
| `PYTHONPATH=scripts/site python3 -m unittest test_build.BuildTests.test_progress_chart_data_file_is_canonical_and_fully_plotted test_build.BuildTests.test_every_scoreboard_panel_has_a_summary_and_the_index_cites_them_all test_build.BuildTests.test_scoreboard_panels_close_and_sit_at_the_top_level test_build.BuildTests.test_scoreboard_is_a_copy_of_the_canonical_repository_file` | 4 passed |
| `python3 scripts/site/test_build.py` | 47 passed; 1 pre-existing error in the lab-browser cover digest check |
| `python3 scripts/site/panel_index.py --check` | 300 panels current |
| `pdfinfo`, `pdftoppm`, visual inspection | Result PDF 3 pages, valid-run PDF 3 pages, wrapper-failure PDF 1 page; all opened and labels/figures checked |
| `git diff --check` | no whitespace error |

The full site test error is
`test_lab_browser_index_is_current_and_self_consistent`:
`scripts/build_lab_browser.py` raises
`ValueError: cover report schema or registry digest is stale`.
`docs/curves/covers.json` records registry SHA-256
`89b580c52ecf14765c61312a80fb233d5ba2758e9db6c147e959bd3d4b18d3a0`;
the unchanged `docs/curves/registry.json` hashes to
`5520bc7e09da329603b75538182ab9c66f756f0347b2b6391cf5456b8daa4925`.
Both values are identical on `origin/main`. The cover record is a
preverified mathematical artifact; changing its digest without rerunning its
native checks would misrepresent its scope. This error is outside the S3
root-batch change. All other 47 site tests passed in the final full-suite
run.
