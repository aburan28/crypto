# Next measured boundary: large-prime independent replay

The current BMSS optimization series does not price or optimize the named
P192/P224 extension-kernel search/replay pipeline. The earlier campaign's
surviving October 12 prefixes now have 294 P192/2341 and 557 P224/1117 sealed
steps, each including the old 50-step prefix. The native integrity audit
reports 244 + 507 = **751 additional maps**, with 1,249 still required.
The walks lack successful child exits and are not complete: P192 is
`IN_PROGRESS`; P224 stopped with `RESOURCE_MONITOR_UNAVAILABLE`.

| Prefix | Construction ms, operational only | Independent replay ms | Summed invocation ms |
|---|---:|---:|---:|
| P192/2341, 294 sealed steps | 1,267,805 | 11,942,474 | 13,218,864 |
| P224/1117, 557 sealed steps | 1,673,303 | 9,280,044 | 10,971,953 |

Source: `isogeny-campaign status` against the frozen 100-map study and
`research/p192_p224_2000_more_large_prime_isogenies_20261010/runs/restart-20261012/walks/`.
These are retained invocation receipts from shared-host, interrupted walks,
not isolated CPU speedup evidence. No duplicate prefix is credited again.

The polynomial replay implementation already has independent Karatsuba
products and reciprocal reduction. In `verification_poly.rs::rem`, the
reciprocal of the reversed modulus is constructed on every reduction.
`verification_kernel.rs` repeatedly calls modular products with the same
kernel modulus. Reusing a prepared reciprocal is an implementation hypothesis;
its share of actual runtime has not yet been measured.

The next experiment must freeze the exact old constructor maps and producer
and verifier sources/binaries. Measure decoding, exact kernel/division-polynomial
checks, codomain/map identity, subgroup and 20 scalar transports separately,
including preparing any reduction context in replay totals. Use the same map
bytes and native Linux isolation on each CPU profile, five A/A pairs and five
interleaved A/B pairs; record all failures and intermediate work.

Before an optimized verifier is accepted, require equality with the retained
independent verifier on each raw map and rejection of mutated kernels,
codomains, rational-map coefficients and subgroup/generator data. Constructor
and verifier arithmetic must remain independent; no constructor success flag
or cached verifier pass may replace algebraic replay. Preserve the old maps,
receipts and acceptance conditions byte for byte.

The requested 2,000 additional maps and multi-CPU whole-process timing remain
open. This protocol is a next experiment, not a replay speedup result.
