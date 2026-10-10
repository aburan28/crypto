# n83 F6 five-summand factored S3 circuit feasibility gate

Registered before implementation and measurement. This follows #1448,
which left 14,382 rows after exhaustive private cubic peeling. The
one-round polynomial matrices are too wide to yield a source-only
constraint. This gate tests a different representation of **the same
four S3 field equations**: a bit circuit retaining multiplication and
squaring factors, Tseitin CNF, and a bounded native Kissat solve. It is
an attempt at an ordinary complete point decomposition, not a claim
that SAT is a finished F6 algorithm.

Use the registered K0 curve
`icv1-f2m83-tm6151469093347-debefd74`, standard dimension-18
factor base (261,447 geometric points; 261,444 subgroup-usable; 130,722
sign-folded columns), five source x coordinates with 18 bits each,
three 83-bit intermediate x coordinates, and four exact S3 links.
Match the source-bit and intermediate-bit layout in
`System512::build`, n=83, b=1, polynomial-basis field table. The
ordinary public T001 subgroup point is
`(355fb5df7a905f16921eb,5900a390f42d290f1bbe)`; take its
verified subgroup preimage under multiplication by four and test all
four rational 4-torsion offsets in order 0,1,2,3. The planted control
uses source indices `[0,2,4,6,8]` and their exact full-group sum.

Build field gates for `S3(x,y,z) = (x+y)^2 z^2 + xyz + (xy)^2 + b`.
Use the same `WideFieldTable::product_bits` and `square_bits` as the
expanded polynomial system. Introduce one Boolean gate output for each
nonconstant AND or XOR, with standard exact Tseitin clauses, simplifying
constant inputs and identical literals. Assert all 332 S3 output bits
zero. In the pinned planted control, add units fixing all 90 source
bits to the five chosen factor-base x coordinates. In the ordinary
runs, leave those source bits free. Before accepting a CNF, evaluate
its circuit on the planted assignment and eight deterministic
assignments and compare each of 332 output bits with the corresponding
expanded `System512` polynomial. Record gate/variable/clause/literal
counts and SHA-256 of each complete DIMACS input. Keep one generated
CNF per arm or a lossless compressed copy and its hash.

Run local `/opt/homebrew/bin/kissat` version 4.0.4 with binary SHA-256
`05d6f3e9c402a1fe8853b0746e384e1b3d1c4a550e255f11daa2461d279aa848`.
First run the pinned planted control (30-second process limit). Then
run each ordinary offset once (120 seconds per process), retaining SAT,
UNSAT, timeout, OOM and nonzero statuses. Use one process at a time,
one solver thread, live 7-GiB RSS sampling/kill and no other local
builds during timed solves. Preserve raw solver stdout/stderr and
exit codes. Parse any satisfying assignment in native Rust, independently
check every expanded polynomial, map each source x to an enumerated
factor-base point, enumerate its 32 sign choices, and require exact
full-group equality to the lifted target. A SAT model without that
full-group witness is algebraic evidence only, not a relation.

Structural success is at least one ordinary full-group-verified
five-summand decomposition under the resource limits. A timeout is
inconclusive for that offset; UNSAT means the exact CNF has no model
for that offset, subject to the encoding checks. If a verified ordinary
witness appears, follow with repeated ordinary queries and complete
F4/F5/F6 same-input measurements before any 2× stage claim; only a
recovered and independently verified DLP on a frozen one-target IC
workload plus paired same-point rho under accepted isolation may support
an end-to-end IC speedup. Until then candidate ID and speedup are
unknown. Timings from this contended Mac are feasibility diagnostics.
