# Static wide-S4 SAT natural-query yield: 6 verified relations in 32

The [fresh-seed, 32-query panel](panel.json) completed exactly once, in trial
order. Every Rust source export was valid. Source-receipted static
CryptoMiniSat reported SAT on six queries and UNSAT on 26; there were no
timeouts, crashes, invalid formula models or nonlifting models. The
[independent post-run audit](RESULT.json) checked **all clauses and XOR rows**
for each returned model, lifted each to the exact geometric factor base, and
group-added its three points to the public query. All six are verified
relations. A separate exhaustive three-sum oracle found **exactly those six**
queries mathematically feasible and all other 26 infeasible. The solver's
UNSAT output has no independently checkable proof trace; the group oracle
supplies the independent mathematical infeasibility check for this small base.

| Fresh ordinary queries | Exact feasible | Verified SAT relations | Feasible missed | Other outcomes |
| ---: | ---: | ---: | ---: | --- |
| 32 | 6 | 6 | 0 | 26 reported UNSAT; 0 timeout/error/nonlifting |

The observed useful-relation yield is **6/32 = 18.75%**, with a descriptive
Wilson 95% interval of **8.89%–35.31%** for a Bernoulli-rate interpretation.
The exact-feasibility rate is also 6/32 on this sample. These intervals are
wide; 32 queries do not settle performance at larger fields or under other
base policies. This is a new query-law seed (`2026092935`) with no overlap
with the previous 104 n17a1 F5 ordinary queries, and the panel was fixed
before any new SAT or exact-group outcome was seen. The earlier six
outcome-selected correctness controls are not counted here.

SAT child-process wall summed to **554.353049 s** across all 32 attempts;
source-export child-process wall summed to **0.869639 s**. These are stage
process sums, not a complete-pipeline elapsed interval. The median SAT
process wall was about **6.014 s** for the six witnesses and **18.836 s** for
the 26 source-UNSAT reports. Per-attempt wall, CPU, peak RSS, stdout/stderr,
source instances and failure status are retained. Stage timings here are
not paired with F5 per-query wall intervals, which the older worker did not
record. No timing ratio or speedup is asserted.

The full [341-file raw archive](evidence.tar.gz) is SHA-256
`e760b906acb441aeefd964f30ff3bb36b0be28a698afab84148350328a56648c`.
The [machine-readable audit](RESULT.json) is SHA-256
`0fcd6e548d32a8d8c5d56b31c0232b12e4edc6bb3e33d8aa495b95002738b4db`.
The [replay test](../../test_static_cms_s4_natural_evidence.py) rechecks the
archive, source/build receipts, all 32 row statuses, every SAT assignment,
all point witnesses and the exact group labels. It also preserves all 26
zero-yield attempts in the denominator.

This is a **PDP stage sample**, not a complete `IC1` candidate run. No
full-rank relation matrix, factor-base logarithms, unseen-target descent,
target recovery, paired rho, one-target online interval or speedup was
measured. The next protocol must integrate this static SAT backend with
relation collection, rank/linear algebra and one-target recovery, then pair
it with F5, the incumbent and rho on fresh public points and the same
resource envelope. The current sample cannot answer that comparison.
