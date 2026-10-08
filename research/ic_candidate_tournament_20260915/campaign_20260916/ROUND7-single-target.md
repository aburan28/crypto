# Round 0007 pre-registration: beat rho on the single-target panel

Written before any measured stage of this round was run.

## Objective

Round 0006 reached the contract's parity rule on complete cold
single-target jobs: 0.8264 of matched rho's instructions with the whole
interval below one, and 0.977–0.985 of its native process time with the
replay interval reaching 1.006. Parity is not a win. This round's declared
aim is to be **strictly below rho in both metrics**: candidate/rho upper
paired 95% limits and every curve-cell ratio below one, in user-space
instructions and in native process wall, on confirmation and on replay.
That is `beats_rho_strict` in the evaluator, and the round is prepared with
`--objective rho`: the incumbent gate only guards against regression
(instruction ratio at most 0.98 with its upper limit below one, native upper
limit below one, no cell more than 10% worse), because the native metric is
dominated by process spawn and shared curve construction that no candidate can
remove, and a 20% margin over the round-0006 winner on the whole job is not the
question being asked. Promotion requires both gates on both final stages.

## Parent, seed, budget

Incumbent: round-0006's promoted `tiny_batch1` — the frozen
`runs/round-0006/source_candidates/tiny/source` with `batch_trials: 1`.
Fresh seed 2026091607; target count 1; pilot profile; 1,800 paired-job
budget; one pinned CPU; 8 GiB cap; 60-second watchdog; Valgrind 3.22.0
`Ir`; native progress recorded. Rho is the shipped per-target
signed-Frobenius solver on the incumbent's executable, as in every round.

## What the native metric is made of

Measured before this round on the round-0006 confirmation fixtures, medians
over 40 to 60 fresh processes pinned to one CPU. A job that only constructs
the curve and target takes 4.0–5.5 ms of wall time, of which about 3.5 ms is
process creation and loading. Rho's own solve adds 0.42–1.28 ms per cell;
the round-0006 winner's own work adds 0.29–0.8 ms (per-phase timers, whose
stderr writes inflate each phase by roughly 10 µs in this VM). Primitive
costs in a warm loop: field inverse 43–66 ns (Euclid) against 165–305 ns
(Fermat), multiply 9 ns, orbit key 18–31 ns, affine addition 106–149 ns.
Per-process cold costs are the surprise: a CPUID instruction traps to the
hypervisor here at about 10 µs each, the standard library's feature cache
costs about 70 µs on first use, libm's first call about 10 µs, and cold code
and data pages cost more than the arithmetic they run. The shared target
hashing initialises the feature cache for both arms; what a candidate can
still avoid is everything else it touches only once.

## Candidates (four arms, 1,488 paired jobs)

1. `incumbent` — as above.
2. `tiny2` — the round-0006 source with:
   - **work-minimising table rows**: the folded pair table builds the `t`
     rows (1 ≤ t ≤ K) that minimise the expected total
     `t·|F| + (K + 3)/(λ·cov(t))`, `λ = |F|(|F|+1)/(2r)`,
     `cov(t) = 1 − ((K−t)/K)²`: entries built plus rests scanned per probe.
     A fixed function of the public base size, column count and order. On
     these cells it is one row (n13, n17, n19) or two (n23), against seven
     or eight before; probes scan a few more rests, which is cheaper than
     the rows they replace;
   - **López–Dahab projective scalar multiplication** for arbitrary points
     and a **doubling table of `G`** for its multiples, with the
     certification of all column logarithms batched behind one field
     inversion per bit (`[ℓ_o]G = [h]C_o` still checked for every column,
     twice per job as before);
   - **carry-less squaring** (one multiplier latency instead of a bit
     spreading chain), the module's **own reduction tables and nibble-table
     normal basis**, a **bit-exact scalar reimplementation of the seeded
     ChaCha12 sampler** (checked draw for draw against `StdRng` in a test;
     the base point set is unchanged and tested equal), no libm calls and
     no hashed containers in the job.
3. `tiny2_rows` — ablation: only the table-row rule on the round-0006 source.
4. `tiny2_arith` — ablation: everything in `tiny2` except the row rule (the
   round-0006 witness-density rule kept).

Development evidence before freezing (not a claim): on the round-0006
confirmation fixtures all 60 `tiny2` outputs pass the frozen checker with
the incumbent's factor-base fingerprints; on one case per cell under
Callgrind the instruction ratio to rho was 0.73–0.77; the per-phase wall
timers put its own work at about 0.19–0.5 ms per job against rho's
0.42–1.28 ms, which projects every cell below one in native time with the
smallest margin on n13a0.

## Boundary, floor, class, honesty

Unit and boundary unchanged. Base support, `m = 3`, no direct relations,
every verification obligation, the worker's phase dumps, fixture
construction, rho branch and general-arithmetic final check are unchanged
from round 0006. The K-instruction floor is unchanged and weak. Class:
engineering. Rho is the shipped implementation; the report will say that its
per-job cost on these cells is mostly fixed setup and cold-process cost, that
a rho specialised the same way has not been measured, and that the CPUID
trap cost is a property of this virtual machine, so native ratios here do not
transfer to bare metal. No arithmetic-complexity, family-wide or
cryptographic-size claim follows. Fresh fixtures from the new seed; every
failure retained; the single-target and 16-target panels stay separate.
