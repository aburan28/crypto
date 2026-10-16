# Mixed-orbit SAT cap arithmetic for the frozen N83 first trial

This is a **source-derived preflight**, not a retained-base solver run. At
source commit `1bd12a72fe315fc86b5ab2df8c3071f6db6d32c5`, the bounded
mixed-orbit worker calls `MixedOrbitDomains::search` with field degree 83,
three summands from 1,182 small orbit representatives, two from 866 residual
representatives, an ordered first small pair, 150,000 maximum variables and
12,000,000 maximum domain clauses. The K values and exact prefix relation
were admitted by the earlier [900-second import gate](../mixed-orbit-import-900s-20261012/RESULT.md).

`frobenius_orbit_domain_size` in
`src/cryptanalysis/binary_semaev_chain_sat.rs` gives, for one summand and `k`
representatives at `n=83` with seven phase bits:

- variables `2k - 1 + 83 + 7 + 2·83·7 = 2k + 1,251`;
- clauses `83k + 3k - 3 + (128 - 83) + 4·83·7 = 86k + 2,366`;
- native XOR rows `(7 + 1)·83 = 664`.

| Domain | Summands | Variables | Domain clauses | Domain XOR rows |
| --- | ---: | ---: | ---: | ---: |
| Small, K=1,182 | 3 | 10,845 | 312,054 | 1,992 |
| Residual, K=866 | 2 | 5,966 | 153,684 | 1,328 |
| **Both** | **5** | **16,811** | **465,738** | **3,320** |

`chain_s3_max_variables(83, 5)` gives **84,743** for the all-affine core.
`ChainedS3Encoding::build` checks that core before construction;
`constrain_frobenius_orbit_domains` checks actual core plus the 16,811 domain
variables; `order_first_pair` checks another `n-1 = 82` variables. Therefore
the all-affine upper total is **101,636**, leaving 48,364 variables below the
150,000 cap. The exact domain-clause count leaves 11,534,262 below its
12-million cap. The worker can examine up to `2^(m-2) = 8` intermediate
identity patterns within one trial, with the global 5,000-conflict cap.

These arithmetic margins rule out the *declared variable and domain-clause
preflight caps* for this frozen 3+2 shape, assuming the imported representative
counts and unchanged source. They do not bound full solver clauses, peak memory,
wall time, model yield, partial-relation yield, graph closure or total
index-calculus runtime. An actual retained-N83 search still requires the
separate 2,400-second authorization in the
[current protocol](../mixed-orbit-search-12600-source-20261012/PROTOCOL.md).
No new experimental seconds are charged by this calculation.

The [machine receipt](calculation.json) freezes the source hashes, integer
inputs and outputs. It was evaluated by direct integer arithmetic from the
named source formulas and independently checked against the two cap
inequalities; no S3 object or worker binary was opened.
