# Preregistered rotated m5/m6 exporter semantics and capacity gate

Status: frozen **before this gate's enumeration**. This is a bounded semantic
reference for a later direct S6/S7 and recursive-S3/S4 exporter. It is not an
S6/S7 polynomial export, a solver benchmark, or an ECDLP result. The parent
[16-target point corpus #767](https://github.com/aburan28/crypto/pull/767) is
merged at `c77767c4a653734f428e110cca29985d721476b2`; its immutable raw
archive has SHA-256
`39f16990213e2c7525c6dae92c131127b12f914bbfd6182264249da72bc9546c`.
`INPUTS.json` pins the exact eight `Q` points and their four `Q+T` counts in
each arm, plus the frozen target/factor file hashes. No target is chosen here.
The source and input SHA-256 values are pinned in `FROZEN.json` before the
first full gate run. This protocol and the input manifest are the preregistration.

## Mathematical question and reference

For `E: y²+xy=x³+1` over the two #767 fields, use their distinct rotated
`V_i=span{beta^(2^i),beta^(2^(m+i))}` and **all rational lifts** in each
`F_i`, with no point at infinity in a factor. The exact projected PDP for a
fixed subgroup point `Q` asks whether some labelled tuple has full sum
`Q+T` for any of `T∈E[4]={O,(0,1),(1,0),(1,1)}`. Keep the four torsion
branches distinct. The direct semantic reference is the set of every
`m`-tuple of two-bit subspace masks for which at least one rational choice of
factor-point signs has that exact full sum. It is built by enumerating all
3,125 or 117,649 labelled rational point tuples, not by trusting an S6/S7
zero or a solver result. It is a **verification oracle**, never an input to a
solver benchmark, because giving a solver its answer table would leak labels.
Compare all 64 exact `Q+T` multiplicities (16 points × four torsion branches)
with #767's archived complete group-law counts. Every point and factor list
must independently match the hash-pinned #767 verifier source and archive.
For every mask tuple in the direct reference, retain the exact rational lift
count and one full rational point-tuple witness; re-add that witness and check
its factor membership, mask decoding, curve equation and full point sum. Every future solver
model must map back to a mask tuple and then to a rational point-sign witness;
a raw algebraic root with no such lift is spurious. A claimed UNSAT must agree
with the complete direct reference, and timeouts/caps remain unknown.

## Recursive-chain and exceptional branches

The prospective generic chain has `m−2` intermediate field abscissae for
prefix sums and `m−1` S3 equations, including the terminal fixed-`Q+T`
relation. For each enumerated rational tuple this gate computes every prefix
by the independently replayed group law. Where the prior prefix and next
prefix are both affine, it verifies
`S3(x(prefix),x(next factor),x(next prefix))=0` using the repository's binary
`S3=(a+b)^2 c^2+ab c+(ab)^2+1`. It classifies ordinary, doubling,
inverse-to-infinity, and identity-prefix additions. A prefix at infinity has
no x-coordinate; its identity or inverse step requires a separate branch,
not an arbitrary fake x. At exact target hits it also checks the terminal S3,
or the explicit identity-prefix equality branch. This establishes **only
necessity on rational witnesses**. It does not prove that all algebraic chain
roots lift to rational points or that a cleared-denominator chain is complete.
A later direct and recursive polynomial exporter must preserve these branches
and prove equal projected solution sets against the direct reference for all
16 fixed point targets before any solver ranking.

The toy layout preflight is fixed before outcomes. Direct x-only summand
variables are `m·d=10` and `12`. A raw S3 chain adds `(m−2)n=39` and `76`
intermediate bits, hence `49` and `88` total **before** sign, infinity,
field-equation, exceptional-denominator or SAT auxiliaries; it has `52` and
`95` binary S3 rows before such auxiliaries. The in-process `F2BoolMono`
`u64` mask accepts n13 chain's 49 bits but rejects n19 chain's 88. Neither
this count nor a direct 10/12-bit layout proves that an expanded S6/S7 ANF is
small enough. Do not run FES on unquadraticized degree-three-or-higher S3
rows, or call a default degree-3 crossbred negative a refutation when input
rows are omitted. An expanded multiword or streaming direct exporter remains
a separate gate.

## Trace and symmetry controls

[Kosters–Yeo, Prop. 4.2](https://arxiv.org/html/1503.08001v3) gives the
ordinary-binary-curve morphism `φ(O)=0`,
`φ(P)=Tr((x(P)+a2)/a1²)` with kernel `2E`. Here `a1=1,a2=0`, so this is
`Tr(x(P))`; the normal beta has trace one and each x mask contributes its
bit parity. Check on **every** rational tuple that XOR of all mask bits equals
`φ(full sum)`. Since `Q∈[4]E`, its trace is zero. For each fixed `T`, the
necessary linear constraint is therefore mask parity `=φ(T)`: 0 for
`O,(0,1)`, 1 for `(1,0),(1,1)` on these odd-degree fields. Verify all 64
frozen `Q+T` branches and count the two parity classes. The all-T projected
question retains both parity classes; this does **not** halve its total tuple
space or prove a 2× PDP speedup. [Kosters–Yeo, Remark 4.8](https://arxiv.org/html/1503.08001v3)
warns that the trace relation is not generally implied by higher Semaev
ideals, because algebraic roots can use nonrational lifts. It is a sound
explicit filter for this rational-point PDP and deliberately changes the raw
algebraic candidate set. A later SAT with/without-trace arm needs separate
frozen encodings and charged construction and solve costs.

An equal-cardinality repeated-`F0` semantic control enumerates the same
number of labelled point tuples on each rung and records projected support
and counts on **these same 16 targets**. It is not a timing reference. The
rotated nonzero x-spaces are pairwise disjoint, so sorting the `x_i` or their
slice-local masks is not a valid symmetry break. The `m!` WDSat DPLL reduction
of [Trimoska–Ionica–Dequen, §4](https://eprint.iacr.org/2019/313.pdf)
assumes interchangeable same-base slots and cannot be multiplied into a
rotated-support estimate. Preserve an explicit swapped-nonzero counterexample
and reject that ordering for rotated solver arms unless a whole-instance,
target- and witness-preserving automorphism is proved. The concrete negative
control is selected **only for structural audit**, outside the 16 frozen solver
targets. For each of two naive order rules (polynomial-basis integer x0≤x1
and slot-local two-bit mask0≤mask1), enumerate all rotated tuples and select
the numerically smallest nonzero projected R outside the frozen 16 R labels
that (a) has no tuple satisfying that order in its entire four-torsion
coset, and (b) has a witness with both first x-coordinates nonzero. If no
such R exists, report that control as censored and make no concrete witness
claim for that order. Otherwise record its full m-point tuple, full sum S,
projected R, subgroup Q, torsion T, and the tuple with its first two points
swapped. Check the swap retains S but puts both nonzero points outside their
new Fi slots; check the sorted projected multiplicity is exactly zero. This
structural Q is never added to the solver corpus. A future matched
rotated/repeated solver comparison must charge work per useful verified
relation and independent rank row, not only compare support fractions.

## Frozen caps, accounting and decision

Run Python ≥3.12 in two sequential cold children, one per arm, with a
600-second wall and 512-MiB peak-RSS acceptance cap **per child**. The runner
keeps each child's JSON, stdout, stderr, exit code, UTC interval, hashes, and
failure output; no retry replaces a failed attempt. Record field and curve
operation counters, CPU, wall, RSS, exact tuple counts, branch counts,
trace-parity counts, and direct mask sets. These are stage correctness costs,
not a common operation-equivalent S or rho ratio. #767's independently
implemented bit-serial/Fermat group law is pinned and reused as this gate's
reference; its original producer used different arithmetic. This gate does
not claim a third independent group law.

**Pass** only if all frozen hashes and 16 target rows match, all 120,774
rotated and equal-count repeated tuples are enumerated within caps, every
rational prefix obeys its ordinary S3 or explicit infinity branch, every
point tuple obeys the trace morphism, all 64 exact coset counts match #767,
and the direct mask sets contain no negative-target witness. Each concrete
symmetry control passes only if its preregistered structural witness exists;
otherwise record it as censored without blocking the 16-target semantic
reference. A single mismatch,
source drift, timeout or RSS overrun is a preserved failed/censored gate,
never a solver UNSAT result. A pass admits a later direct/chain *exporter*
comparison, not SAT/F4/F5/WDSat performance claims. The next exporter PR
must freeze exact polynomial/circuit bytes and hashes, degree/monomial and
auxiliary counts, explicit denominator and infinity cases, model lift,
per-engine variable/degree/memory preflight, 15-second/2-GiB pilot caps from
[#763](../ROTATED_M5_M6_SOLVER_ADMISSION_20260925.md), and all 16 point
labels before outcomes. A solver ranking follows only if both encodings are
semantically exact on the toy corpus. Before any actual ECC2K-130 target
experiment, verify the public challenge field/model against #762's n131
basis (and construct an isomorphism only if they differ), import its literal
source `P,Q`, verify curve and subgroup order, and archive coordinates,
source-model hashes and import receipts. This toy gate uses no challenge
coordinates.
