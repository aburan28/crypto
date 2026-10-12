# Frozen native build before point generation

The protocol was committed and opened as draft PR #1353 before any pilot or
held-out point was generated. This build completed from clean source commit
`e37c50e9dbbc88148a12be648735894c37ba0b1a` with:

```sh
cargo build --offline --locked --release \
  --example koblitz_orbit_dlp_fast_online \
  --example koblitz_rho_fixture \
  --example koblitz_one_target_replay
```

The command exited 0 in 1m35s. The toolchain was Homebrew `cargo 1.93.1`
and `rustc 1.93.1 (01f6ddf75 2026-02-11)`, on `arm64` macOS 26.6.
The existing library warnings were unused mutability, an unused variable,
unconstructed AVX lane variants and unread batch-scratch fields; none was a
build error. The host did not permit `sysctl -n hw.model`, so this receipt
does not assert host-level CPU isolation or identify a physical CPU model.

`BUILD_SHA256SUMS` pins `Cargo.toml`, `Cargo.lock`, the three executable
example sources and their locally built release binaries. The Git commit
pins the complete transitive source tree. Rebuild with the exact command and
compare hashes before treating a later binary as this experiment's producer.
The release binaries and Cargo target directory are local reproducibility
artifacts, not source-controlled inputs; the hashes remain in Git.
