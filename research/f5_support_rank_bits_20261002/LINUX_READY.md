# Width-6 physical Linux gate prepared; result pending

The workflow `.github/workflows/f5-support-rank6-call.yml` implements the
frozen direct selective-echelon versus cut-19/width-6 complete-call gate.
It rebuilds on an x86-64 GitHub Linux runner, checks AVX2 and BMI2,
reserves/pins one or two CPUs, keeps up to three attempts per seed, runs
four frozen seeds, applies the existing complete-call analyzer with an
explicit width-6 assertion, and uploads qualified and failed receipts.
The unchanged cut-19 workflow still defaults to width 8. The benchmark
and analyzer now record and check the selected candidate cap.

Local validation completed:

- Release GF(2) tests: 9 passed; release F5 tests: 14 passed.
- Native paired runner and analyzer built on Apple ARM64. `actionlint`
  accepted the new workflow; `shellcheck` accepted the updated shell
  orchestrator.
- The four Linux workflow examples passed an offline x86-64 Linux
  `cargo check` using the installed Rust 1.98 rustup compiler. This is
  a compilation check, **not** physical x86-64 execution. An initial
  attempt selected Homebrew Rust 1.93 and failed because that compiler
  lacked the target's standard library; the accepted rerun used the
  installed rustup target and succeeded.
- A single Apple ARM64 smoke comparison, `frozen`, completed all 22
  processes with exact output and cut-19 original rows. Width 6 used
  24,233,029 reduction word XORs against selective echelon's
  100,213,183. Its paired complete-call median was 2.696×, but the
  exact bootstrap lower bound was only 1.373×, and A/A ranged
  0.602–1.522×. This is nonpromoting diagnostic timing.

The full smoke receipt is `LOCAL_DIRECT_frozen.json.gz`; raw and
compressed SHA-256 values are in `LOCAL_DIRECT_RAW_SHA256.txt` and
`LOCAL_DIRECT_ARCHIVE_SHA256.txt`. The physical Linux four-seed
one-thread >2.00× gate and two-thread/smaller-case controls remain
**unmeasured**. A further 2× wall-time result is therefore unclaimed.
