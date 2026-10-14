# Infrastructure and rejected attempts

1. Homebrew Cargo rejected `cargo +1.98`; the corrected command is
   `rustup run 1.98 cargo ...`. No benchmark output was produced.
2. The first unmodified release build in
   `/private/tmp/isogeny-performance-20261014-base` failed with LLVM
   `IO failure on output stream: No space left on device` and exit 101.
   A queued second build was stopped with SIGTERM before producing a
   baseline, to prevent further system-disk exhaustion. Neither is a
   mathematical result or a timing sample.
3. Automatic approval review initially rejected removal of an old root-gate
   target directory because its contents had not been inspected. Read-only
   inspection found only Cargo-generated release artifacts; the passing test
   log is preserved separately in the committed earlier study. `cargo clean
   --release --target-dir <that cache>` then succeeded through approval review,
   removing 1,025 files / 293.6 MiB of compiler output. Run evidence stayed in
   `research/p192_p224_2000_more_large_prime_isogenies_20261010/evidence/`.
4. Builds were moved, with automatic approval review, to the task-specific
   `/Volumes/SSD990-2/crypto/isogeny-performance-builds-20261014/` on the second
   SSD (811 GiB available at inspection). Compiler temporary files use that
   directory as well. The source worktree stays in the authorized repository.
5. Local machine inspection showed several unrelated CPU-heavy jobs. macOS
   cannot supply the Linux core-reservation/PSI contract. Local native wall
   measurements are consequently not promoted; operation counts and exact
   outputs can still be audited while shared jobs run.

No processor-specific optimization has been rejected or accepted yet.
Optimization rungs below will keep regressions alongside improvements.
