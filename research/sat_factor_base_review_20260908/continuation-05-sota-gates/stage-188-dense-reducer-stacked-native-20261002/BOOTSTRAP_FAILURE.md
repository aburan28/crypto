# Stage 188 bootstrap failure

The first attempt to meter the exact candidate build invoked
`target/release/examples/koblitz_f4_stage188` after only
`cargo test --example koblitz_f4_stage188` had run. Cargo had produced the
hashed test executable, not the ordinary runnable example path. The shell
therefore returned exit 127 with:

```text
zsh:1: no such file or directory: target/release/examples/koblitz_f4_stage188
```

No child experiment or build process started, and no process receipt was
available. A one-time unmetered bootstrap `cargo build --release --locked
--example koblitz_f4_stage188` produced the native meter. The subsequent exact
candidate build, exact tests, screen, confirmation, verifier corrections,
composition control, and final verification are metered. Because this
bootstrap is unmetered, complete campaign cost remains `null`; it is not
silently estimated or omitted from that claim.
