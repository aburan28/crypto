# N83 finalist validation record

All commands below used branch `codex/n83-k3072-current-source-20261017`
with Rust source unchanged from finalist producer commit
`416489e09c3d8885b6c12f262935e6ba29373460`. The final local release
library gate was:

```sh
CARGO_TARGET_DIR=/private/tmp/n83-k3072-current-source-20261017/target \
  cargo test --release --lib -- --test-threads 1
```

It finished with **4,162 passed, 0 failed, 117 ignored** in 1,538.12
seconds. Its compressed [complete log](local_lib_serial.log.gz) has
SHA-256 `400c69a64a4e87683ad4602a135d789fcb4ded900ad0514510664d6641bfac95`.
The two executed examples' release test suites passed: 21/21 for
`koblitz_orbit_dlp_fast` and 4/4 for `n83_rho_optimization_eval`.

Earlier local gate attempts remain classified by their recorded status.
The initial sandboxed parallel library run had 4,149 passes, 13 failures
and 117 ignored. Nine `pollard_collab` tests could not create loopback
sockets (`PermissionDenied`); four `jv_cover` tests hit their time caps or
coverage assertions under the high parallel load. An unsandboxed
four-thread rerun had 4,161 passes, one `jv_cover` Nagao timeout, and
117 ignored. The Nagao case then passed alone in 4.60 seconds and also
passed in the complete serial suite above. These are environment and
parallel-load test observations, not arithmetic verdicts for the N83
candidate.

The completed [cloud run](https://github.com/aburan28/crypto/actions/runs/39005111289)
had two successful job conclusions. Each provider ZIP was downloaded
and hashed. `aggregate.py` independently checked every extracted member
and replayed the producer receipts locally; both local verifier exits were
zero and the 128 IC, 128 rho and four repeat records yielded
`PASS_INDEPENDENT_REPLAY`. The reconstructed
[deduplicated raw archive](raw_dedup.tar.gz) passed `verify_archive.py`
for all 3,802 file paths. The canonical materializer produced eight
candidate manifests, 132 IC1 rows and 128 paired rho rows.

The five-arm follow-up's `test_protocol.py` passed its fixture,
100-events-per-host schedule, candidate identity and selector controls.
The workflow passed `actionlint`; actual preflight for both halves
resolved to the declared K1120/root11 fallback with 25 IC, 25 rho and
50 calibration events per host. `python3 scripts/site/test_build.py`
passed all 48 site tests after regenerating the panel index.
`python3 scripts/build_ic_leaderboard.py --check` and
`python3 scripts/build_lab_browser.py --check` both passed for the
generated leaderboard and browser outputs. The result PDF rendered as two A4
landscape pages and both pages were visually inspected.
