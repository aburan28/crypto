# Frozen n37 descendant-native bounded residual search

Status: **preregistered; no b01 outcome has been measured under this protocol**.
This tests whether a fixed set of target-independent shifts closes the
107/1024 direct-query support gap measured on the disjoint b00 block. It is
an explicitly batched correctness and cost question, separate from the
repository's primary one-target IC-versus-rho comparison. No batch result is
an attack-speed or ECC2K-130 feasibility claim.

## Hypothesis and frozen inputs

The descendant-native signed 3+3 oracle with the same 42-column factor base
will recover all 1024 previously unseen public Q on the orbit-disjoint b01
block using a direct query followed, only on a complete direct miss, by no
more than 16 fixed shifts. The falsification condition is any target without
a full-point-verified logarithm at the cap, any incorrect support/miss
decision, or rank below 42 after 256 target-blind relation probes. Whether
the method is competitive with rho is a later, separately measured gate.

- Prerequisite is the n37 native six-sum PR #1254, at corrected head
  `357f643a1945bf018f6eb15afaeec173b97f0639`. Its b00 result is
  diagnostic evidence, never an input to offset selection or the b01 solver.
- Use the exact `NATIVE42.json` base, SHA-256
  `bd4bd8af982bcc65234ae1fae1ad52eb0ed41a1c203bf4d2556fd8e87f788a4c`,
  and the degree-73 archive, SHA-256
  `eb4773d556886672e8b80735fc486c16b73765fd5069ff81133100bfc3eafa90`.
  Reconstruct and validate the field-basis and isogeny bridge, all 42 base
  points and source pullbacks, and the transported subgroup generator.
- Use `disjoint_cold_v2_20261001/FROZEN.json`, SHA-256
  `da958a3f1117dd1b88703fed2055c5c9f64255516e1e920f35d193bf53c6a88d`.
  The public b01 point file is
  `research/notes/ecc2k130/disjoint_cold_v2_20261001/fixtures/n37_L1024_b01.points.jsonl`,
  SHA-256 `78553fdff5ae66521d3a1052978962258e78d48c92285df6e1df962034b43717`.
  It contains 1024 public points from corpus
  `compact-disjoint-cold-v2-n37-L1024-b01-20261001`, seed
  `2026100110101`. The solver must not read any known-answer scalar.
  Independent replay may read the separate b01 fixture only after the
  producer has written its answers; fixture SHA-256 is
  `3453995c0423e4911ad4a6afa7cfe50ec1d6b0e6c67a771bda590f1d6ed680df`.
- The subgroup has order `r=230603167`; its full curve order is
  `137439487532` and cofactor is 596. Freeze the source and archive field
  representations exactly as in PR #1254.

## Frozen algorithm

Rebuild the complete signed at-most-three half table from the 42 leaf base
points in PR #1254's enumeration order. Query with its exact full-point
3+3 MITM oracle. Re-run the same target-blind relation probe stream:
SplitMix64 state `0x6e33376d365f3031`, scalar `1 + output % (r-1)`,
discard repeats, cap 256, stop at rank 42. Solve the independent 42-row
system and verify all 42 base logarithms by scalar multiplication before
querying b01 targets. Charge the entire setup to this cold run.

Generate one **global, Q-independent** ordered list of 16 distinct nonzero
shift scalars from SplitMix64 state `0x6e33375f72657331`, using
`1 + output % (r-1)` and discarding repeated scalars. Precompute all 16
`[a_j]G_leaf` points, charge every scalar multiplication, and preserve the
ordered list in the raw result. For each public target in input order:

1. Transport Q and query it directly. If it has a verified six-summand
   witness, derive `log Q` from the solved base logs and stop for that Q.
2. On a complete direct miss, query `T(Q)+[a_j]G_leaf` for `j=0..15` in
   order, stopping at the first verified witness. Every failed query is a
   recorded complete miss. Derive `log Q = log(T(Q)+[a_j]G_leaf)-a_j (mod r)`.
3. Verify each reported answer against both `[log Q]G_source=Q` and
   `[log Q]G_leaf=T(Q)`. If all 16 shifts miss, record an unresolved target
   with no inferred scalar. Never use the b01 fixture to fill a miss.

The direct and residual queries must use the same exact oracle and table.
Record each queried leaf point, shift, witness or complete miss, complete
per-query counts, and any verification error. The independent replay must
reconstruct source pullbacks; form `Q+[a_j]G_source` for each recorded shift;
derive exact at-most-six membership using iterative source-curve reachable
sets, not the producer's MITM lookup; verify every decision and witness;
solve relation logs independently; verify every recovered scalar on the
source curve; and compare to the b01 fixture only after these checks.

## Accounting, reference, and decision

The cold accounting interval starts before reading inputs and ends after
writing complete results. Include input hashes, bridge/isogeny work, base
selection and transport, all 16 shift scalar multiplications, half-table
construction, relation probes and failed lookups, dense solve, direct and
residual failed/successful queries, group additions, scalar verification,
and output. Preserve exclusive phase times, process wall/user/system time,
operation counts, peak table size, rank trajectory, each target's attempt
count, maximum attempts, and unresolved count. A replay run is a separate
verification cost, not part of attack cost. Keep single-process wall times
descriptive and leave `S`, speedup, and matched rho fields null.

The first decision is **complete verified batch or bounded-residual no-go**
on b01. A complete batch permits, but does not prove, a performance gate:
pair cold native residual runs with a same-host, same-Q signed-Frobenius
batched-rho arm using `koblitz_rho_batch_ks_v3` at
`37 0 signed_frobenius 1024 2026100110101`,
`KIC_RHO_POINT_INPUT` set to the b01 public point file,
`KIC_RHO_BATCH_CORPUS=compact-disjoint-cold-v2-n37-L1024-b01-20261001`,
normal-basis canonicalization, default 32 parallel walks and fixed
distinguished-point configuration. Before any runtime comparison, freeze
the exact rho environment/configuration, run an A/A noise control, then
interleave at least five uncontended matched cold A/B rounds. Report all
failures and both methods' complete costs in the same calibrated unit with
paired uncertainty; do not promote a stage time or batch throughput to the
primary one-target speedup. A separate one-target, same-Q online IC/rho
experiment under the repository's five-phase measurement contract is the
required primary attack-speed gate.

Commit the producer source, exact build lockfile, commands, raw decisions,
timing, independent replay receipt, analysis, and scoreboard decision in a
follow-on PR, whether the hypothesis succeeds or fails. Do not infer n41,
n53, n131, or Certicom-scale behavior from this n37 result.
