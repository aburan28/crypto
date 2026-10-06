# P-256 Montgomery-affine pair table, round 36: protocol

Date registered: 2026-10-06

## Question

Round 35 reduced 384-MiB projective lookup overhead to `0.680998` complete
P-256 additions per segment, still above the `0.334411` budget.  Can a
directly usable 64-byte Montgomery-affine entry remove one third of table
traffic and avoid canonical-coordinate parsing while preserving exact group
replay and a rho-fair arithmetic denominator?

## Frozen dependency and scope

Import round 35 by exact SHA-256
`90704e913a0407fe773ad7550cd5b68c5ead54077ff485978bd51cdf21bcfc89`.
Keep cutoff 219, all 17 variable columns, eight signed-pair lookups, eight
setup additions, the target correction, mean capacity
`224.361249307511`, corrected base ratio `0.998569034150286`, and allowed
overhead `0.334410508423` additions per segment unchanged.

The complete affine table has `34,562,148,612 * 64 = 2,211,977,511,168`
bytes, below `2^50`.  This round changes only entry representation and the
matched addition primitive; it does not change the support or relation model.

## Registered representations and fair denominators

1. `projective96-prefetch16`: reproduce round 35's directly addable
   projective control.
2. `mont-affine64-complete`: store `(x,y)` as two Montgomery field elements,
   synthesize `z=1`, and use the existing complete projective addition.
3. `mont-affine64-mixed`: use a specialization of the same complete RCB
   formula with the affine operand's `z=1` substitutions performed before
   compilation.  Count and report the field multiplications removed.

For each affine candidate, the denominator is eight additions using the same
affine representation and the same complete or mixed formula with operands in
an L1-resident array.  The projective control retains its matched projective
denominator.  Thus an arithmetic shortcut available to the selector is also
available to rho and cannot be counted as selector-only progress.

Sweep prefetch distances `{1,2,4,8,16,32,64}` for both affine candidates.
Prefetch both 32-byte coordinates' cache line(s) and charge every instruction.
Use `2^15`, `2^20` and `2^22` entries.  The `2^22` affine decision row is
256 MiB; it contains the same number of pair entries and the same index stream
as round 35's 384-MiB projective row.

## Exactness and timing

- Generate the same 4,096 hash-selected nonzero P-256 points as round 35.
- On at least 4,096 segments, require complete-affine and mixed-affine output
  to equal both the projective control and an independent reference.
- Exercise identity, doubling, inverse and distinct-point cases separately;
  no incomplete mixed formula is admissible.
- Require exactly eight reads and eight additions per segment and zero group,
  curve, range, false-positive or false-negative failures.
- Hash pool, indices and canonical outputs.

Run one warmup and nine deterministically rotated release repetitions over
65,536 segments.  Record wall and process CPU time, checksums, host, compiler,
entry/table bytes and every raw row.  Use sorted ranks 1, 4 and 7 as the
interval.  Convert each candidate's table overhead with its matched direct
denominator:

```text
extra matched-add equivalents = 8*(table_time/direct_time - 1), clamped at 0
projected/rho = 0.998569034150286
              + 0.964336477130181*extra/(224.361249307511+1).
```

Select the smallest median projected ratio on the `2^22` row; smaller
prefetch distance breaks exact ties.  The timing gate requires both median and
rank-7 upper endpoint no greater than rho.

## Gates and stop condition

Advance only if the selected affine row has zero exactness failures, complete
storage below `2^50`, and both registered timing endpoints at or below rho.
Do not infer 2.212-TB deployment behavior from 256 MiB; a passing row remains
a stage candidate until a defensible full-table model exists.

All end-to-end gates remain: proved P-256 usable-relation probability,
structured residual degree at most 5, relation collection below `2^120`, cost
per usable row below `2^103`, complete measured time, and non-genericity.
Attempt no full-depth relation unless every gate passes.  If the affine rows
fail, publish the smallest measured gap and reject these representations
without claiming that local degree or cache-only timing establishes parity.
