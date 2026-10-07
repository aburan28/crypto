# n83 F6 basis-aware S3 circuit gate

Registered before implementing or measuring this follow-on to #1452. The
native-XOR factored circuit passed a fully pinned planted control but found no
ordinary model in four 120-second limits. This experiment tests whether
keeping the dimension-18 source coordinates in their original basis during
field multiplication makes the exact SAT instance materially smaller or
enables a verified ordinary relation. It is a point-decomposition stage gate,
not a complete F6 or IC speed comparison.

Freeze the curve `icv1-f2m83-tm6151469093347-debefd74`, standard dimension-18
factor base, five source x coordinates, three free 83-bit intermediate x
coordinates, the four exact S3 links and the public T001 subgroup point from
#1452. Use its planted `[0,2,4,6,8]` full-group control and the same four
cofactor-four preimages obtained by rational 4-torsion offsets 0, 1, 2, 3.
Preserve #1452's circuit source, raw results and original emitter as the
reference. The exact frozen reference binary SHA-256 is
`2261bd409f6754cb23242a2dfcda20e8a1f2a5d3efbb2873e2032fa6b848f619`.

The candidate uses the characteristic-two identity

`S3(x,y,z) = (xy + xz + yz)^2 + xyz + b`.

Represent every source x as its 18 input bits paired with its actual field
basis elements. Form each product of a source and another operand from those
18 terms directly, preserving every field coefficient exactly. In each link,
reuse `xz` or `yz` to form `xyz` with a source operand. Dense intermediate
operands may retain the old 83-bit multiplication. Keep the original 339
input-bit layout and assert all 332 output bits zero. Require the candidate
circuit to agree with the expanded `System512` on the planted assignment and
the same eight deterministic arbitrary assignments before emitting any input.
Record AND, XOR, variable, mixed-constraint, byte and SHA-256 counts for both
representations on matched inputs; these are encoding diagnostics, not speed.

Use the same CryptoMiniSat 5.14.7 executable and SHA-256 from #1452, one
thread, the same runner's 7-GiB live RSS kill gate, and sequential runs. First
run the fully pinned planted control with a 30-second limit. If it has an
independently verified full-group witness, run the source-only pinned planted
control for 30 seconds, then each ordinary offset for 120 seconds. Retain all
SAT, UNSAT, timeout, OOM and process-error rows and raw output. Any SAT model
must pass independent expanded-polynomial and full-group verification, and
ordinary relations require subgroup-usable source points. Freeze source and
input hashes and preserve compressed inputs with exact uncompressed hashes.

Structural success requires a verified usable ordinary relation within its
limit. Secondary success is at least a 20% drop in AND gates and mixed
constraints on each ordinary input, with all equivalence checks passing;
this alone does not establish solver speed or relation yield. If all ordinary
offsets time out, record a bounded negative for this exact encoding and
limits, with no claim of UNSAT. Any 2x F6 stage claim requires repeated
complete same-input decompositions against F4 and F5; an IC speedup requires
a verified one-target online run and same-point rho under accepted isolation.
Until those gates, candidate ID and online speedup remain unknown. All CPU
wall times on this contended host are exploratory.
