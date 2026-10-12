# Decision: K16 retains its n37 implementation lead, rho remains cheaper

The preregistered whole-solve instruction census retained the 16-column
factor base as an **implementation-cost lead over the fixed 42-column base**.
Its eight same-point cold `ecbench` solves used 0.36712 times K42's
Callgrind instructions, with a target-resampled 95% interval
[0.36481, 0.37006]. All 32 profiles recovered the archived public-point
logarithms. The duplicate K42 arm agreed to within 0.0014% on every point,
well inside the preregistered 2% A/A gate.

| One-target arm | Exact candidate or method ID | Actual usable base points | Folded columns | Verified profiles | Mean solve Ir | Mean Ir / √r |
|:--|:--|--:|--:|--:|--:|--:|
| Strong signed-Frobenius rho | `ECM1h1d9961ee4601` | — | — | 8/8 | 9,351,187 | 615.79 |
| IC K16 | `IC1N37Ckb0fb1184PDP3mitmfrobeniuscountedRCsampleLAgaussTDpdpISO0h44f5af6dc772` | 1,184 | 16 | 8/8 | 41,711,016 | 2,746.74 |
| IC K42 | `IC1N37Ckb0fb3108PDP3mitmfrobeniuscountedRCsampleLAgaussTDpdpISO0he8d240ba37b3` | 3,108 | 42 | 8/8 | 113,616,995 | 7,481.88 |
| Identical K42 control | same K42 candidate | 3,108 | 42 | 8/8 | 113,617,055 | 7,481.88 |

The table's only cost unit is **Callgrind simulated user-space instructions
inside `methods::solve`**. `Ir / √r` is labelled
`crypto.S.callgrind_ir`, with `r = 230603167`; it is not group-addition
equivalents and has no comparison to the generic-group floor. The eight
target-paired ratios of sums are K16/K42 = **0.36712** [0.36481, 0.37006],
K16/rho = **4.46051** [3.75421, 5.45372], and K42/K42-control =
**0.99999947** [0.99999525, 1.00000444]. The brackets are fixed-seed,
20,000-sample target-block bootstrap percentile intervals; the per-target
ranges and exact counts are in [DECISION.json](DECISION.json). K16's 63.3%
reduction from K42 in this unit supports the frozen column decision. It is
still 4.46 times rho's instructions. There is **no IC/rho crossover or
admitted wall-time speedup** in this result.

The profiler was Valgrind 3.22.0 on an Ubuntu 24.04 x86-64 hosted VM
reporting an AMD EPYC 9V45 CPU. The release `ecbench` binary SHA-256 was
`7e90fcf5520355ffd1d1352e0ce172127a803908947c73db8471372c0313a5a5`.
The selected CPU feature path inside Valgrind is not independently recorded,
so these counts describe this simulated execution, not isolated native CPU
time. The boundary begins after input, workload and method reconstruction;
the solve includes factor-base and rank preparation, target descent, scalar
replay, and the adapter's repeated factor-base inventory. The runner's
external independent audit is outside the counted interval. This is a
conservative implementation census and leaves the primary one-target online
wall metric unset. No n41, n53 or n131 transfer is inferred.

The [frozen protocol](PROTOCOL.md) and [32 exact jobs](JOBS.json) were pushed
at `963aca96` before [CI run 37186853441](https://github.com/aburan28/crypto/actions/runs/37186853441).
The [raw CI archive](RAW-CALLGRIND.zip) is 2,895,137 bytes, SHA-256
`809f19a0b784e7ecf08f74c372a16d8a6afe4aff817eb1d2d8168ce676901420`.
It contains every child input/output, Callgrind part, per-file SHA-256 receipt,
and Valgrind log. Its compact [census](CENSUS.json) SHA-256 is
`e09fefcdaff8373c0f25f79b2b2f78a35cc3092caf2e134914c50481a1a08415`;
the [host and build provenance](PROVENANCE.txt) SHA-256 is
`9316d3774461085216e8b687c20d1f61bbd77045baf6d324cf78cdd3bf7d5ed8`.
The original GitHub artifact is ID `11297686371` and expires on
2027-01-02; the committed ZIP preserves the raw evidence beyond that date.
The [native analyzer](../../examples/ecbench_callgrind_analyze.rs) re-hashes
every raw file, re-parses all Callgrind parts, checks the child method,
workload and recovered scalar against the frozen archived records, and then
recomputes [DECISION.json](DECISION.json). The archived n37 session's 320
measured runs had already passed an independent Linux replay; no profile
failed or timed out.

This settles the K16-vs-K42 **instruction** question. K8 was statistically
close to K16 in the earlier incomplete GAE count and was not profiled here,
so the next base-size gate is a preregistered K8/K16 complete-cost comparison
on untouched points. The larger research gate remains a physically isolated,
same-point n41/n53 one-target online wall comparison with all native work and
preparation accounted for separately. The factor-base/isogeny track should
compare equal usable sizes before inferring an advantage from a descendant
presentation.
