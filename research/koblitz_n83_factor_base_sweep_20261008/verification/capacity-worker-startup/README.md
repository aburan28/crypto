# Guarded Linux worker startup check

The native example was cross-built as an `x86_64-unknown-linux-musl` executable
with Rust 1.98.0 and Zig 0.16.0, using the package's `--no-default-features`
mode. `sat-capacity-linux-build.log` records the build; the executable's
SHA-256 is `a67c516b6914fe25399e5e8d759163280695d46d272c9835e0320456e50a10c5`.
The normal package default still enables Redis native TLS. The isolated
capacity worker does not use Redis and omits that TLS feature so cross-build
does not require an unrelated target OpenSSL sysroot.

The supervisor launched this executable with a 256 MiB cgroup memory cap,
zero swap, no network, the same pinned image ID as the synthetic guard, and a
5-second wall cap. The deliberately empty
`empty-panel/manifest.json` caused a JSON `EOF while parsing a value` error
after the worker's cgroup check. The outer status is correctly
`PRODUCER_FAILURE_worker_exit`, exit code 1. This check proves guarded
startup and malformed-input classification only. It did not import retained
points, construct an S4 model, search, or measure an index-calculus runtime.

`config.json`, `outer.json`, `stdout.log`, and `stderr.log` retain the exact
invocation and outcome. The earlier cross-build failures with the Homebrew
compiler's missing target and the native-TLS OpenSSL dependency remain in
`../sat-capacity-linux-homebrew-target-failure.log` and
`../sat-capacity-linux-openssl-failure.log`.
