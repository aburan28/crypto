# Local socket permission during library verification

The first mandatory full release run used the restricted default command
sandbox. It completed in 247.79 seconds with 3,987 passes, nine failures and
112 ignored tests, exit 101. All nine failures are mock-network tests in
`pollard_collab::cairn` and `pollard_collab::net`: `TcpListener::bind("127.0.0.1:0")`
was denied with `Operation not permitted`. Raw failures are retained in
[`records/verification/library.log`](records/verification/library.log).

No code, test or ignore flag was changed. The same compiled test binary's
`cryptanalysis::pollard_collab` checks were run with sandbox network access
enabled, giving 37 passes, zero failures, zero ignored, exit zero. The same
mandatory full Cargo command is repeated with that network permission in
[`records/verification-network/`](records/verification-network/).
It passes 3,996 tests, zero failures and 112 ignored, exit zero, in 238.05
seconds. Every catalog test name is accounted for; unfinished/extra lists
are empty. The restricted failure remains recorded separately.

Reproduction of `verify.sh` requires permission to create loopback sockets.
Each full run retains its separate 1,200-second budget, commands and actual
exit code; the unsuccessful restricted run is not relabeled as a pass.
This is a test-environment correction, not a change to experiment workloads,
acceptance thresholds or the implementation.
