# First validation attempt lacks a terminal receipt

The current-base `cargo test --release --lib --locked --offline` attempt wrote
`root-release.log` and reached compilation. After execution handles disappeared,
resuming handle 99942 returned Unknown process id. At 2026-10-16 06:29 UTC,
recorded launcher/cargo PIDs 69166 and 69168 were absent. No exit-code receipt
exists. This attempt is incomplete infrastructure evidence, not a passing test
suite or a mathematical result. Its log is retained without modification.

A fresh retry uses the same unchanged root source and separate log/exit/PID
files. The process is detached from interactive tool handles so its final receipt
can survive a goal continuation. The existing build cache remains an optimization;
its contents do not substitute for the retry's actual terminal test result.

The ordinary `nohup` retry also lost launcher PID 15442 without an exit-code
receipt. Its separate log and supervisor bytes remain unchanged. A generic
native launcher now creates a new OS session with `setsid` before exec. A
15-second control survived the launching tool call, was observed live, and wrote
its terminal completion marker. The following retry uses this verified process
separation. This changes process custody only; scientific commands and their
built-in resource/correctness gates stay as frozen.
