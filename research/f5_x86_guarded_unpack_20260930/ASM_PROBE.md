# Representative x86 release code-generation probe

The probe was compiled on Apple ARM64 with the installed Rust 1.98.0 x86-64 Linux cross-target. It is a standalone representation of the checked and guarded term-copy loops, **not** the full F5 crate. Rebuild it with:

```sh
rustup run stable rustc --crate-type lib --edition=2021 --target x86_64-unknown-linux-gnu -C opt-level=3 -C target-cpu=x86-64 --emit asm -o asm_probe_x86.s asm_probe.rs
```

The exact [probe source](asm_probe.rs) has SHA-256 `fb37a99d48db2e016b3024a69d4dbbb180ed64c9470bf53171ba2a07362f25f2`; the committed [assembly](asm_probe_x86.s) has SHA-256 `d92bf6c616485fac98603d710ba710a73b232af15e7a8a744cf09ecc8342ab4c` from `rustc 1.98.0 (88d9e12ae 2026-08-18)`. In `checked_decode`, the per-term loop around `rep bsfq` contains `cmpq %rcx, %rsi` and `jae` to `panic_bounds_check`. In the validated path of `guarded_decode`, the per-term loop around the second `rep bsfq` reaches the column load without a bounds comparison. Both paths still enumerate and copy one value per set bit.

This supports an x86-specific instruction-count hypothesis only. Compiler optimization of the full crate, host CPU latency and the complete F5 call must be measured under the frozen protocol before claiming any speedup.
