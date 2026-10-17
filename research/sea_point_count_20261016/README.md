# SEA and point counting: measurements and formal checks

The BSGS point-count workloads improved by paired median ratios of 1.98652,
2.79019 and 3.24310 at 24, 32 and 40 bits. The measured SEA and full P-192
API effects lie within their A/A noise envelopes; no SEA speedup is established.

The formal work also led to a real CM curve-generation correction: sampled
annihilation could accept the wrong twist. Generation now validates CM norm
data and requires exhaustive counting or successful Hasse certificates.

- [Results and requirement status](REPORT.md) / [PDF](REPORT.pdf)
- [Frozen protocol](PROTOCOL.md) and [measurement accounting](MEASUREMENT.md)
- [Run report](reports/run-20261016T0205Z/REPORT.md) / [PDF](reports/run-20261016T0205Z/report.pdf)
- [Raw measurements](evidence/run-20261016T0205Z/runs.jsonl), [native summary](evidence/run-20261016T0205Z/summary.json), and [reference certificates](evidence/run-20261016T0205Z/reference-certificates.json)
- [Verification map and remaining obligations](VERIFICATION.md)
- [Validation commands and logs](VALIDATION.md)
- [TLA+ model](formal/Recombination.tla), [Lean proofs](lean/SeaProofs.lean), and [final checker logs](formal/logs/final/)
- [CM regression before](evidence/pre-fix-cm-regression.jsonl) and [after](evidence/post-fix-cm-regression.jsonl)

Replay formal checks with `bash formal/verify.sh`. For timings, install Rust
1.99.1, then run `bash replay-benchmark.sh`; allow up to 60 minutes. The script
uses frozen source worktrees for both arms, rather than measuring the later
CM constructor correction. See the protocol for limits and cache preparation.
Native benchmark/analysis sources are in `native/src/bin/`; rendering recorded
results with `python3 render_reports.py` performs no experiment or analysis.

This is partial formal verification of the search argument. The concrete Rust
arithmetic, SEA congruences and analytic CM construction are not fully proved.
