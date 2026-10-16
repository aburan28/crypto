# Post-merge validation and consumer controls

Recorded on 2026-10-15 from source
`0191f5b924765a0f6d7882f13787b7391413c685`, including main through
`1ad1d2836cbaeb49a39a1e203a0729c413ade1f7`. Subsequent changes in this
publication commit affect workflow triggers and evidence/documentation only.
The arithmetic implementation is unchanged. These test intervals and controls
are excluded from performance ratios and from the separate 2,000-map novelty
count. Earlier frozen archives are unchanged.

| Receipt | Command boundary | Observed result |
|---|---|---|
| `root-post-merge.log` | Full release library suite, two test threads, inherited Rayon configuration | 4,094 passed, one deadline failure, 117 existing ignored; shell exit 101 |
| `nagao-post-merge-exact.log` | Same compiled binary, exact failing test alone, `RAYON_NUM_THREADS=1`, one test thread | One passed, zero failed; shell exit 0 |
| `root-post-merge-rayon1.log` | Full release library suite, two test threads, one Rayon worker | Interrupted with SIGINT at the long concurrent fixtures; shell exit 101; no full-suite verdict |
| `root-post-merge-serial.log` | Full release library suite, one test thread, one Rayon worker | 4,095 passed, zero failed, 117 existing ignored; terminal PASS present; original tool-session exit receipt was no longer retrievable on resumption |
| `post-merge-algorithms.log` | Standalone release library, binaries, integration tests and examples | All executed suites passed; existing exploratory ignored test retained |

The root command used the recorded Homebrew rustc 1.93.1 and durable target:

```sh
TMPDIR=/Volumes/SSD990-2/crypto/isogeny-performance-builds-20261014/tmp \
RAYON_NUM_THREADS=1 cargo test --release --lib --locked --offline \
  --target-dir /Volumes/SSD990-2/crypto/isogeny-performance-builds-20261014/root-gate \
  -- --test-threads=1
```

The Nagao failure reports `timed_out: true` at trial 1. The test and its solver
are unchanged from the original baseline. The successful serial run retains
every assertion and test; no timeout or mathematical acceptance condition was
relaxed. Attribution to interference is a hypothesis supported by the isolated
pass and the concurrent failures, rather than a measured isolation receipt.

The dedicated `isogeny` consumer controls retain their complete raw output,
resource receipts and independent replay under `p192-13-control/` and
`p224-1117-control/`, plus stdout/stderr. The embedded build revision is
`6e90a71a72b7b7c3522bb0307ad1d25299bc98ce`. P192/13 enumerates both rational
eigenlines; P224/1117 constructs one and reports PARTIAL coverage. Independent
replay passes kernel, codomain, exact rational map, subgroup and twenty public
scalar transports for each map. These are correctness controls, without a
qualified CPU timing comparison.

The current coordinator has a different project ID and task inventory; old
T-2844 returns 404. T-15 reserves the resumed publication paths through a
registered work shell. This is a new coordination record after the state
change, with the previous failure history preserved in `ATTEMPTS.md`.

The visually inspected [report](REPORT.pdf) has editable [Typst](REPORT.typ)
and [control-diagram SVG](CONTROLS.svg) sources. Raw receipts are packaged in
`../POST_MERGE_VALIDATION.tar.gz`; restore their original paths with:

```sh
tar -xzf research/isogeny_performance_20261014/POST_MERGE_VALIDATION.tar.gz \
  -C research/isogeny_performance_20261014
```

`FILES.json` applies to this extracted directory. The original raw files also
remain unchanged in the local worktree; packaging does not delete them.
