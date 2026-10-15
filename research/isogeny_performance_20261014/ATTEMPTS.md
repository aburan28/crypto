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

11. A root attempt completed with 4,077 passes and six failures because two
    archived runtime fixtures had not been materialized. The tracked chain-set
    constants and SOURCE42 file were loaded without changing tests. A subsequent
    run stopped while writing test output to the full main SSD. Both logs remain.
    The final unfiltered run wrote to the second SSD and passed 4,083 tests,
    zero failures, 117 existing ignored, in a 1,347.49-second test interval.
12. The main SSD also filled before rung 3 comparison receipts could be written.
    Raw constructor results survived; the native comparisons were repeated in
    fresh output files after space recovered. Git rejected a proposed cross-volume
    worktree move with `Cross-device link`; the source paths stayed unchanged.
    Live test output was directed to the durable second-SSD evidence path.
13. Conductor exhausted the implementation task's four attempts and marked it
    failed after the final root validation succeeded. T-2844 covers the distinct
    publication/CPU-evidence review phase; it reserves evidence/report/index paths
    and does not repeat the completed arithmetic implementation work.

14. The post-merge root run failed the unchanged Nagao six-sum test after its
    existing deadline expired. Its exact rerun passed. A concurrent full retry
    remained stalled at long fixtures and was interrupted, preserving its raw
    output. The final unfiltered serial run passed 4,095 tests, zero failures,
    with 117 existing ignored. See `publication-20261015/` and `VALIDATION.md`.
15. On resumption on 2026-10-15 the coordinator's project ID and task inventory
    had changed, and old task T-2844 returned 404. Clear scopes were claimed as
    T-15 in the current coordinator. A registered interactive work shell now
    maintains the real publication/validation claim. Stale MCP credentials
    returned 401; the configured CLI supplied the current authenticated access.

16. The source SSD filled during a kernel-adapter edit. The file retained its
    previous complete bytes; a second-SSD copy was made before the corrected
    edit. Headroom subsequently recovered without deleting any research data
    or the earlier campaign executable. The first adapter compile also found
    a private helper path and an unimported point-access trait; the adapter now
    uses the public curve helper and an explicit affine-point match. Neither
    attempt produced an accepted benchmark.

17. GitHub CPU wave 38766895864 retained one ARM admission refusal: runner
    housekeeping used 0.29 CPU-seconds during a two-second settling window
    with a 0.20 limit; PSI remained below its limit. Nineteen later commands
    were admitted, but the required full series is censored. The next wave
    uses five settling seconds with the same CPU-rate and PSI thresholds.
18. The source SSD again refused new directories while reporting about 31 MiB
    available. A first audit-archive attempt reported extended-attribute
    temporary-file errors; the fresh no-xattr archive uses second-SSD TMPDIR
    and passes byte-for-byte extraction comparison. No raw records were removed.
    A shared-clone attempt unexpectedly began packing the shallow repository;
    it was stopped before checkout. A lightweight continuation checkout on the
    build SSD reads the existing Git objects through an alternate, while all
    new objects and evidence stay on that SSD. The old checkout and running
    tests remain intact. Root library bytes and all three modified native
    files were compared before continuing; no runtime result is inferred from
    these infrastructure attempts.
