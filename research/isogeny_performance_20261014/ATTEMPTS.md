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

6. A launch rejected a revision argument with the wrong length before writing
   a run. A subsequent 192-bit fixture launch panicked because it used four
   Montgomery limbs; the backend requires three. Those stderr files and the
   header-only interrupted receipt remain in `evidence/baseline/`. The corrected
   fresh run uses the minimum limb width and has all 96 cases plus its PASS exit.
7. Root and standalone test builds initially lacked compile-time fixtures in
   the sparse checkout. The failed logs are retained. Only the referenced
   tracked fixtures were materialized before fresh validation runs.
8. The local actionlint version does not recognize the current GitHub
   `macos-15-intel` label. The label is listed in GitHub's runner reference;
   no other workflow issue was reported by that invocation. This remains a
   tooling-version limitation until validation with a current label catalogue.
9. The report helper initially linked a library with the wrong compiler.
   Cargo was invoked through `rustup run 1.98 cargo`, but Cargo's compiler
   fingerprint and the retained `.rlib` show Homebrew rustc **1.93.1** on PATH.
   The helper uses that compiler, and the actual compiler receipt is retained.
   Formatting and Clippy invocations used their rustup toolchain drivers.
10. An attempt to transfer the Conductor claim used `task release`, which
   cancelled T-2829; this server would not reopen a cancelled task. T-2834
   transparently replaces it. Its lease is now attached to this task's real
   registered root-validation process and maintained by that process.

The first hoisting rung has 96 exact output matches. Its 224-bit/degree-101
linear-algebra case has 273,659,435 versus 5,639,435 multiplication calls
(48.526 quotient). Wall time remains unqualified. Later rungs retain any
increase in multiplication calls beside reduced inversion counts.
