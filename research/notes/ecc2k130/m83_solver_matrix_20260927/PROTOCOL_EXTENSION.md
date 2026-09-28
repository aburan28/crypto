# Preregistered extension: native-XOR SAT and payload-size sweep

Commit this extension before running any new cases. It supplements, and does
not modify or overwrite, the six-variable run in `PROTOCOL.md`.

## Hypotheses and fixed inputs

1. CryptoMiniSat 5.16.0 using shared Tseitin AND variables and native XOR
   clauses may beat Z3 on the *same* fixed-phase `m=83` Boolean systems.
   Run it on all original 12 fixed-phase inputs (two seeds × six inputs);
   keep the first pass and its source/target hashes intact. Enumeration must
   produce the complete algebraic model set, then the same independent group
   verification. Solver seed 13, one thread, five-second aggregate check
   limit and the original process budgets apply.
2. At the same `E0/GF(2^83)` modulus, subgroup, normal-coordinate seeds
   `260938` and `260939`, and fixed phases `(0,1,2)`, increase payload bits
   from 2 to 3 and 4. At *each* size and seed, use the first lexicographically
   selected proper planted triple and the first subgroup target drawn with
   `random.Random(seed*1000+830+payload_bits)`; pass identical target and ANF
   digest to FES, Z3 SAT, native-XOR CryptoMiniSat SAT, Boolean Macaulay/F4
   and Boolean F5B. Make no choice after observing solver outcome. Keep s=2
   rows in the original run as the reference. If the factor base cannot form
   a planted target, record setup_failed for that cell.

## Reference, success, limits

FES scans `2^(3s)` assignments: 512 for s=3 and 4096 for s=4. Every claimed
solution must match FES's complete algebraic root set and must survive lifting,
the subgroup order check and group addition to the target. Record entire
receipts, backend versions, input/source hashes, peak memory, all work phases,
counts (F4 rows, F5 pairs, SAT decisions where available), and absent or
timed-out cells. Stop each case after 25 seconds for equation construction,
5 seconds for solver setup/check, 70 seconds total and 1 GiB address space.
FES completion is required to adjudicate correctness; if FES itself times
out, leave that cell's oracle verdict unresolved. A positive solver timeout
does not imply no decomposition. Do not replace failed inputs or increase a
limit after the fact.

For native XOR, one positive variable per payload bit, one auxiliary for each
distinct degree-2+ monomial, three CNF clauses for an AND when degree two
(generalized to k+1 clauses for a k-input AND), and one XOR clause per
Boolean equation; block models on the original payload bits. The empty
monomial is represented as the XOR clause parity, not a new free variable.
Log the numbers of clauses, XOR clauses, auxiliary variables, and decisions
when exposed by the backend.

These extra sizes are fixed-phase, planted-control solver-stage probes; even
if an engine solves rapidly, `83^3` Frobenius phase combinations, natural
relation yield, factor-base rank, full linear algebra, final DLP and matched
rho all remain outside these timings. Thus no ECC2K-130 improvement follows
from a fast cell. Where full-cost matched units are unavailable, leave
cross-solver operation ratios and end-to-end `S` null. Keep this extension's
results in a separate immutable results directory with their own manifest.
