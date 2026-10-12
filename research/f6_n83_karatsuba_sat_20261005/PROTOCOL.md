# n83 F6 recursive carryless S3 circuit gate

Registered before code or measurement. The basis-aware native-XOR circuit
in #1453 reduced ordinary mixed constraints by about 26% but all four
ordinary lifts timed out. This follow-on tests whether the two remaining
products of fully free 83-bit intermediate coordinates can be represented
with fewer nonlinear gates. It is a bounded encoding and decomposition
feasibility test, not a full F6 or one-target IC comparison.

Freeze the same curve `icv1-f2m83-tm6151469093347-debefd74`, exact
dimension-18 base, public ordinary T001, four 4-torsion preimages, planted
source indices `[0,2,4,6,8]`, five summands, 339 input bits and 332 S3
output equations as #1453. Its release reference binary SHA-256 is
`ccf7fdedf814b2464b9f920079d04d8bfb1648be3feddc18fb8478793d122e01`.

Preserve the #1453 emitter as the reference. In the candidate, replace
only the `xz` dense-by-dense product in S3 links two and three (zero-based
links 1 and 2) with recursive polynomial-basis carryless multiplication.
Pad both 83-bit operands with zero bits to 128. At each recursion level,
compute low and high halves and the XOR-sum half, using three half-size
products; reconstruct the exact 255-coefficient unreduced polynomial and
reduce each monomial with `WideFieldTable::product_bits` under the pinned
field modulus. Simplify constant and identical Boolean literals as in the
reference. No other group or polynomial equations change.

Before any SAT run, require each candidate circuit to match every output
of expanded `System512` on the planted assignment and the same eight
deterministic arbitrary assignments. Record each candidate and reference
input's AND, XOR, variable, mixed-constraint, byte and SHA-256 counts.
Run exactly the same CryptoMiniSat 5.14.7 executable and SHA-256 from
#1453, one thread, the same 7-GiB live RSS kill gate, 30-second fully
pinned planted control, 30-second source-only planted control and four
120-second ordinary offsets in order. Stop if the fully pinned control
fails independent expanded-polynomial and full-group verification. Retain
all raw solver statuses, including timeouts, errors and OOMs. Any ordinary
SAT model must pass the existing independent expanded-equation and exact
full-group witness checks, including subgroup usability.

Primary structural success is at least one verified usable ordinary
relation under its registered limit. Secondary encoding success requires
at least 25% fewer AND gates than #1453 on every ordinary input, no
increase in mixed constraints, and all equation checks passing. Gate
counts alone do not prove solver speed. If all ordinary lifts time out,
preserve them as inconclusive, never UNSAT. A 2x F6 stage claim still
requires repeated complete same-input decompositions against F4 and F5;
one-target IC speed requires a verified recovered scalar and paired
same-point rho under accepted isolation. Keep candidate ID and IC
speedup unset until those complete gates.
