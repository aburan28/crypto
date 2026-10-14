# Release test follow-up

The full release library command `cargo test --offline --locked --release --lib` ran once on the shared host. Its immutable [raw log](full-lib-test-final.log.gz) records 4,068 passed, 13 failed and 117 ignored tests in 823.65 seconds. Four `cryptanalysis::jv_cover` failures occurred during high parallel host load; nine `cryptanalysis::pollard_collab` failures reported sandbox `Operation not permitted` when opening loopback sockets. These are test-environment failures in that run, not evidence about the relation collector.

The exact affected groups were rerun without changing source:

| Command | Environment | Result | Raw log |
| --- | --- | --- | --- |
| `cargo test --offline --locked --release --lib cryptanalysis::jv_cover::tests:: -- --test-threads=1` | serial on the same host | 14 passed, 0 failed, 1 ignored; 54.43 s | [jv-cover-serial.log.gz](jv-cover-serial.log.gz) |
| `cargo test --offline --locked --release --lib cryptanalysis::pollard_collab:: -- --test-threads=1` | loopback sockets allowed | 37 passed, 0 failed; 0.73 s | [pollard-collab-loopback.log.gz](pollard-collab-loopback.log.gz) |

The separate [equal-prefix test](fixed-prefix-fresh.log.gz) passed with the frozen fresh inputs. Its [post-merge replay](fixed-prefix-postmerge-replay.json) passed and reproduced the full parsed result object exactly under the newly built test binary. The focused [lane-stream test](final-lane-test.log.gz) passed widths 4, 8 and 16 and checked exact points, coefficients and restarts. The full library suite has not been rerun in one combined process under both corrected conditions; its original failing receipt remains preserved.

After merging the existing PR head and main, the [prime-lane unit test](prime-lanes-postmerge-preformat.json) passed again. Formatting-only source edits happened while that release compilation was in progress, so its receipt states the invocation commit and the source timing explicitly. The separate [final-tree unit test](prime-lanes-final-tree.json) passed at source commit `431d0b3167027092443705e21d363489b2637a22`. The [release ecbench build and plan receipt](postmerge-build-check.json) records binary SHA-256 `ab5764921ad616ae8177ca57c125ec8cd7ba09896f1e0644a6999d1fb33b7159`; its eight arms, eight workloads and 256 execution count match the sealed lane session's plan. The build includes main through `5ca51e183`. A later main merge arrived during validation and touches a file reserved by another active task, so it was not folded into this test tree.
