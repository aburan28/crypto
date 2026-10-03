# Fixed cut-19 matrix-F5 complete-call candidate

## Decision

The fixed cut-19 implementation passed the frozen four-seed local
exactness and ≤35% counted-work screen against selective echelon. The
requested further-2× **isolated physical Linux x86-64 complete-call
result is still unmeasured**. These Apple ARM64 timings are nonpromoting:
each exact five-pair bootstrap lower bound is below 2× and A/A ratios
are broad. One smaller-case timing control fell below its A/A minimum.
This is a solver-stage result, not an IC online or matched-rho speedup.

The source uses cut 19 directly, without an experimental environment
selector. The earlier cut sensitivity and cut-20/cut-19 paired studies
are preserved on `codex/f5-support-cut-screen-20261002` at commits
`2f43c977e` and `57740ed4d`. Those studies showed cut 19 certifies all
four frozen matrices and reduces its exact certificate XORs by about
16.8% against cut 20. This branch builds on the independently committed
cut-20 implementation `9a5c835af` and keeps its exact fallback.

## Matched local screen

One release binary ran both forms in separate processes on the same
frozen input, one Rayon thread. Each seed had two warmups, five
reference/reference A/A pairs, and five alternating
selective-echelon/cut-19 pairs. All 22 processes per seed completed and
matched exact rank, canonical row space, F5 criterion, built/pruned
counts, and column count. The cut-19 certificate succeeded on every
primary call and returned original rows. All six smaller cases matched
every nontiming output and route field.

| Seed XOR | Complete-call selective echelon / cut 19 median | Exact bootstrap 95% lower | A/A ratio range | Reduction word XORs: reference → candidate | Candidate work / reference |
| --- | ---: | ---: | ---: | ---: | ---: |
| `00000000` | 1.515× | 0.779× | 0.791–1.436× | 100,213,183 → 28,129,589 | 28.07% |
| `badc0de1` | 3.252× | 1.463× | 0.464–2.084× | 100,342,346 → 28,198,391 | 28.10% |
| `5eed2026` | 3.283× | 0.707× | 0.263–2.486× | 100,270,257 → 28,114,847 | 28.04% |
| `f5c02a28` | 2.271× | 0.360× | 0.403–6.986× | 100,365,756 → 28,134,474 | 28.03% |

The local work gate passes on every seed, but the complete-call timing
does not prove 2×. The `badc0de1` smaller case `f5_n24_m24_d3`
measured a 0.540× paired ratio against a 0.856× A/A minimum; that
case never took the support-cut route and kept identical output. Its
timing signal is preserved rather than removed.

## Verification and receipts

- Release `cryptanalysis::matrix_f5_f2::tests`: 14 passed.
- `actionlint` passed for `.github/workflows/f5-support-cut19-call.yml`,
  `shellcheck` passed for this directory's `run_linux.sh`, and the four
  workflow examples passed an offline Linux x86-64 `cargo check` with
  the installed rustup compiler. That is a cross-compilation check,
  not execution on physical x86-64 hardware.
- The native analyzer was smoke-tested against the four Apple receipts
  with synthetic clean-isolation records. It rejected them for the
  expected Apple-host, bootstrap, and smaller-case reasons while
  recognizing the fixed cut-19 route. This is a negative gate test,
  not a Linux measurement.
- Four native paired receipts, 22 calls each, are in `local_*.json.gz`.
  `RAW_SHA256.txt` hashes their original JSON bytes and
  `RECEIPTS.sha256` hashes the compressed archives. All four archives
  were extracted and matched their raw hashes. Every receipt has
  source/binary hashes, stdout/stderr, status, output fields and phase
  timings. The benchmark runner source SHA-256 is
  `af111bc46c1e26f37944eabcd6395cfbb5f572eecc9c7bf4051987cb8ed0a1da`;
  the exact matrix-F5 source SHA-256 is
  `bd9ef483e825dffcd7204132199a82e58bca9f4c762ab134c2f12ae75e18374e`.

The next gate compiles and replays this source on isolated physical
Linux x86-64. Each one-thread seed must have complete-call median and
exact bootstrap lower bound both above 2.00× against selective
echelon, with smaller-case and separate two-thread nonregression
controls. Failures, timeouts, and contended attempts remain explicit.
