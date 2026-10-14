# Release test follow-up

The full release library command `cargo test --offline --locked --release --lib` ran once on the shared host. Its immutable [raw log](full-lib-test-final.log.gz) records 4,068 passed, 13 failed and 117 ignored tests in 823.65 seconds. Four `cryptanalysis::jv_cover` failures occurred during high parallel host load; nine `cryptanalysis::pollard_collab` failures reported sandbox `Operation not permitted` when opening loopback sockets. These are test-environment failures in that run, not evidence about the relation collector.

The exact affected groups were rerun without changing source:

| Command | Environment | Result | Raw log |
| --- | --- | --- | --- |
| `cargo test --offline --locked --release --lib cryptanalysis::jv_cover::tests:: -- --test-threads=1` | serial on the same host | 14 passed, 0 failed, 1 ignored; 54.43 s | [jv-cover-serial.log.gz](jv-cover-serial.log.gz) |
| `cargo test --offline --locked --release --lib cryptanalysis::pollard_collab:: -- --test-threads=1` | loopback sockets allowed | 37 passed, 0 failed; 0.73 s | [pollard-collab-loopback.log.gz](pollard-collab-loopback.log.gz) |

The separate [equal-prefix test](fixed-prefix-fresh.log.gz) passed with the frozen fresh inputs. The focused [lane-stream test](final-lane-test.log.gz) passed widths 4, 8 and 16 and checked exact points, coefficients and restarts. The full library suite has not been rerun in one combined process under both corrected conditions; its original failing receipt remains preserved.
