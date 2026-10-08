# P-256 serpentine sign-orbit streaming, round 291: protocol

Date frozen: 2026-10-06

## Requested scope and evidence plan

| continuing requirement | round-291 evidence |
|:--|:--|
| pursue parity with rho | test a table-free global ordering whose registered optimistic stage target is below rho |
| keep all 17 columns variable | retain the complete cutoff-219 unsigned-tuple population and all `2^17` sign assignments per tuple |
| remove the routing obstruction | generate group states directly in tuple/sign/path order; materialize no pair requests and perform no radix routing |
| do not discard branches | visit every sign and every path state for each executed tuple; account for tuple selection separately |
| exact replay | compare exhaustive toy state multisets with direct enumeration and replay registered P-256 transitions in the native group |
| preserve end-to-end gates | leave coverage, same-family degree, relation collection, sparse linear algebra and non-genericity unset unless directly established |

## Frozen identity and dependencies

- curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- measured family: round-33 known-log low-delta/scalar selector;
- comparison factor base: `FB1h2f8621cda105`, not composed with this
  selector;
- columns `B=131458`, arity 17, rare edges `R=6935`, cutoff 219;
- signs per unsigned tuple `2^17=131072`;
- round-33 result SHA-256:
  `9931bccbd9e2821f4498f465ce65f83efb400e7fb3740a239bdb40d3275ffbd8`;
- round-31 coefficient source, transitively bound by round 33, SHA-256:
  `9b9ee16c5a24e868d1b9304f08aa80ca39ac23b1b4586ef730534e739172c435`;
- rounds-43--290 manifest SHA-256:
  `41cdec6ad5265464cc512b2f026f6fa77b333aa209f163e7e77ab78ba665b876`;
- local-oracle ratio `0.964336477130181` rho;
- exact retained-tuple mean capacity `224.36124930751055`;
- round-33 pair-table setup: eight additions and 264 random bytes per sign
  segment;
- Round 108 routing kernel: `0.44960544597285523` native-addition
  equivalents per start, excluding request generation.

The binary must verify both dependency hashes and their curve, family,
capacity, cutoff, digest and non-promotion fields before emitting evidence.

## Candidate construction

For one accepted unsigned tuple, retain its exact round-33 lowest-index
common-edge path of capacity `D`.  Enumerate all sign masks by binary-reflected
Gray code.  Traverse the path forward for even Gray ordinals and backward for
odd ordinals.

The endpoint of one traversal is therefore the start point of the next.  A
Gray transition flips exactly one signed column at its current endpoint.
Precompute the doubled point for each of the 17 endpoint columns and update the
group sum with one addition.  The next path traversal then returns to the
opposite endpoint.  This is the serpentine invariant.

Every `(tuple, sign, path-position)` state from the round-33 model remains
present exactly once.  Tuple boundaries initialize a fresh 17-term sum and
the endpoint doubles.  No pair table, pair request, radix record, discarded
sign or probabilistic branch is credited.

## Registered accounting

For a tuple of capacity `D` and `L=2^17` sign segments, charge:

```text
path additions             L*D
sign-boundary additions    L-1
right-colour corrections   L
initial 17-term sum        16 additions
endpoint doubles           17 native doublings
states emitted             L*(D+1)
```

Report native additions and doublings separately.  For the conservative
stage boundary, charge every doubling as two additions.  Thus

```text
A(D) = (L*D + (L-1) + L + 16 + 2*17) / (L*(D+1))
stage ratio = 0.964336477130181 * A(D).
```

At the registered mean this predicts approximately `0.968616` rho.  It is a
stage hypothesis, not a speedup.  Record emitted states, exact operation
counts, bytes read/written, peak state, and any tuple-selection work.

## Exact controls

1. Exhaust every distinct tuple, every sign and every legal path state on the
   existing round-33 toy groups `(p,B,R,arity)=(257,10,3,3)` and
   `(65537,20,3,4)`.  Direct natural-mask enumeration is the reference.
   Require identical state multisets, multiplicities and digests, with zero
   transition, endpoint, false-positive or false-negative failures.
2. Reconstruct the complete native P-256 coefficient cycle and require the
   frozen digest
   `980917981827d813e60484abb0655e8bd527b0beb540ea146d2974ff17303a53`.
3. Select accepted P-256 tuples using the frozen round-33 hash/rejection rule.
   Execute sign depths `8, 10, 12, 14, 16, 17` on one frozen primary tuple,
   with independent holdout tuples at depths `8, 10, 12`.  Replay every
   reported transition in coefficient arithmetic and the native P-256 group.
   At every sign boundary compare the incremental state with a direct 17-term
   reconstruction at the registered endpoint.  Apply one-unit mutations as
   negative controls.
4. Fit counted-operation, state-emission, memory and measured-time exponents
   over at least four depths.  Timing uses one thread under
   `tools/isolated_bench.py`, with A/A native-addition controls and seven
   rotating baseline/candidate repetitions.  Wall time is secondary.

If the full depth-17 native replay exceeds 20 minutes or 16 GiB, preserve the
partial artifact and mark that cell incomplete; do not transfer exactness from
a smaller depth.

## Boundary table and gates

The result must place in one addition-equivalent table:

- Pollard rho reference;
- round-33 materialized-pair stage;
- Round 108 routing stage;
- direct table-free natural-sign enumeration;
- serpentine Gray streaming, counted and measured.

Promotion requires all of:

- zero false positives and false negatives on every complete cell;
- exact equality with direct enumeration and native group replay;
- counted and measured complete selector at or below rho in the applicable
  tier;
- proved P-256 achieved coverage or usable-relation probability for the same
  correlated stream;
- structured residual degree at most 5 for this measured family;
- complete 138,031-row collection below `2^120` field/group-equivalent
  operations and below `2^103` per usable row, including duplicate/rank
  allowance and sparse linear algebra;
- peak materialized storage below `2^50` bytes;
- a demonstrably non-generic end-to-end algorithm.

Do not attempt a full-depth unplanted relation unless every gate passes.  A
sub-rho operation count with missing coverage, degree or collection evidence
is retained as a stage diagnostic and does not establish parity.
