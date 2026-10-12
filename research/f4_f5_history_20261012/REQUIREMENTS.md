# F4/F5 plot and PDF requirements

Requested: separate F4 and F5 plot groups, graphs and PDFs that retain improvements over recorded revisions and allow comparison on fixed workloads.

| Deliverable | Acceptance condition |
| --- | --- |
| Matrix-step history | Separate F4/F5 series for the existing shared n=23, m=2, degree-4 input; instruction and runtime graphs |
| Solver history | Separate F4/F5 series for the existing shared n=17, m=2 input; instruction and runtime graphs |
| Complete-call allocation history | Separate F4/F5 series from the frozen native baseline/candidate audit; cumulative requested bytes and allocation requests |
| Criterion allocation history | Its own group, since criterion construction is only part of F5 |
| PDFs | One inspected PDF per group, plus a combined report; editable sources and graph data retained |
| Continued tracking | Regenerate groups, graphs and PDFs from recorded snapshots in the existing history workflow |
| Evidence | Preserve the seven imported raw snapshots, source hashes and allocation references; no invented points, intervals or speedups |
| Validity | Break lines across missing data, fingerprint changes, harness changes and host/compiler changes; label historical CI timings as diagnostics |
| Delivery | Commit and extend PR #1677; preserve its existing changes and stop after publication |

Historical snapshots were produced by the existing Python-orchestrated performance pipeline. This task prepares plot data in Rust and uses Python only for figure layout and PDF rendering. It runs no new performance comparison and makes no runtime-speedup claim. The allocation series and historical main-branch series have different workloads and remain separate.

Conductor is unavailable in this environment. The attempted scope check is retained in `records/conductor.log`. Work continues within the authorized feature branch.
