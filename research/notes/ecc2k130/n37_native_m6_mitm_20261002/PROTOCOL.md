# Frozen n37 descendant-native six-sum experiment

Status: **preregistered; no outcome has been measured under this protocol**.
This is one arm of the ECC2K-130 factor-base comparison, not a conclusion
about ECC2K-130 or a claim of a rho crossover.

## Question and boundary

Can the fixed 42-column, descendant-native factor base from the degree-73
transport bridge provide verified six-summand relations of full rank and
recover the frozen n37/L1024 public logarithms? A five-summand control is
mandatory. Its exact formal-support ceiling is
`V(42,5)/230603167 = 37129037/230603167 = 16.10084%` for a uniformly
sampled subgroup target. A finite 1024-target sample can exceed that
population ceiling by sampling fluctuation; report its confidence interval
and never treat the ceiling as a deterministic sample-count limit. For six
summands, `V(42,6)/r = 2.28882687` is only a counting capacity, not a yield
prediction. The reference for attack cost is a cold, same-source/Q,
same-host signed-Frobenius *batched rho* run. The previous source-policy
n37/L1024/off cold CPU ratio was 3.031 to its matched rho, but that historical
number cannot be used as the new native arm's denominator.

## Frozen inputs and selection

- Code starts from merged bridge PR #1253, merge commit
  `d532b3369b43f8dd39b42f5cc0e82c92a8c1a5d1`.
- `NATIVE42.json` has SHA-256
  `bd4bd8af982bcc65234ae1fae1ad52eb0ed41a1c203bf4d2556fd8e87f788a4c`.
  Its 42 leaf points are the first distinct signed classes after cofactor
  596 projection of ascending leaf abscissae. No Q selects or changes them.
- `FROZEN.json` has SHA-256
  `da958a3f1117dd1b88703fed2055c5c9f64255516e1e920f35d193bf53c6a88d`;
  public n37/L1024 block-00 `points.jsonl` has SHA-256
  `187ec04fe50326bbb2f17dadf76056841f04af37a8abb8ff7f502fcd531711ad`.
  The experiment and solver may read only this point file, never its separate
  fixture-scalar file. Independent replay may compare recovered scalars with
  the fixture after the solver has written its answer.
- The archived degree-73 isogeny input has SHA-256
  `eb4773d556886672e8b80735fc486c16b73765fd5069ff81133100bfc3eafa90`.
  The field-basis generator image is 10156182909. The control and leaf both
  have 137439487532 rational points, prime subgroup order `r=230603167`, and
  cofactor 596.

## Exact oracle and relation rule

Construct the complete signed half-sum table from the 42 fixed leaf points:
the identity and every multiset of one, two, or three points from the 84
ordered signed points (`+B_0,-B_0,+B_1,-B_1,...`). Duplicate group sums may
share a key; retain the first multiset in enumeration order and report both
raw entries and distinct sums. For each target, scan distinct half sums in
first-insertion order, in fixed 4096-point batches, and look up `Q - S` in
the table. Stop on the first full-point match. Re-add the at-most-six signed
factors independently before accepting a witness. A complete miss means all
distinct half sums were tested. The five-summand control uses the identical
right-hand table and scans only left-hand sums with at most two factors.
Signed cancellation and repeated factors are allowed; the emitted row is
the *net* coefficient vector and its `L1` norm must be at most the arm's
arity.

Relation probes use `[a_i]G_leaf` for deterministic, target-blind nonzero
scalars from SplitMix64 seed `0x6e33376d365f3031`, reduced as
`1 + output % (r-1)`; discard repeated scalars. The first 256 distinct
probes are the hard cap. Only independently verified witnesses enter a
mod-`r` rank tracker. Stop collecting once rank 42 is reached. Solve the
42-by-42 system, then verify every base logarithm by `[log B_i]G_leaf=B_i`.
Transport each of the 1024 public Q once, run the six-summand oracle, derive
its log from the verified base logs, and verify `[log Q]G_control=Q` for every
answer. Preserve each miss, dependent row, invalid witness, timeout, and
incomplete result. Never reconstruct Q from a fixture scalar.

## Accounting and decision

Cold cost begins before reading inputs and includes checksum validation,
field-basis construction, archived isogeny reconstruction, base/point
validation and transport, half-sum construction and de-duplication, every
failed and successful lookup, relation generation/verification, rank and
linear algebra, every scalar verification, and outputs. Record phase wall
times, total wall time, logical group additions, batch calls, scalar
multiplications, peak table entries, relation probes, misses, rank trajectory,
and exact input/source hashes. A same-host native rho comparison must include
its own full setup, all walks, collisions, and verification and use these
same 1024 Q; without it, leave end-to-end speedup and `S` unset. Do not
promote the earlier 3.031 ratio or a stage-only time to a speed claim.

The diagnostic passes its correctness gate only if rank 42 is reached within
256 distinct relation probes, every accepted relation and all 42 base logs
verify as full points, and all reported target logs verify. The fixed-batch
success gate is **1024/1024** recovered Q; any miss means the bare fixed-base
six-sum method fails that batch and needs an explicitly costed residual
strategy. Compare the five-summand sample yield to its 16.10084% uniform-target
population ceiling with binomial uncertainty; a higher observed fraction is
not by itself an implementation failure.
Regardless of success or failure, commit raw counts, witnesses sufficient
for independent replay, costs, verification receipt, and a decision in the
follow-on PR. The four-policy comparison, n41/n53, rho crossover, and n131
transfer remain separate gates; no result here closes them by implication.
