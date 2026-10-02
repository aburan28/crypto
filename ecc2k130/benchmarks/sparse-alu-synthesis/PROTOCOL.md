# Sparse-basis ALU linear-map synthesis: frozen static screen

Source parent: `b3d0095f6fb74c180c71780ead8ab692b4831cc1`.

This bounded native screen asks whether a materially different ALU circuit can
replace the rejected three-bit shared-table gathers for the sparse
`z^131+z^8+z^3+z^2+1` representation.  It does not compile or run CUDA, build
a sparse walk, search, collect distinguished points, solve a collision, or
recover a scalar.

## Boundaries

The measured fused sigma reference is **15.436677 B complete scalar
updates/s**.  The engineering objective is **26 B/s** on one RTX PRO 6000.
Using the most favourable recorded 188 SMs at 2.430 GHz and 64 ALU
lanes/SM-clock defines a deliberately simple source-operation screen: each
charged source word operation consumes one ALU lane-cycle, with no credit for
compiler lowering, instruction mix, overlap, or latency hiding.  That model
gives complete-walk screening allowances of:

```
reference: 188 * 2.430e9 * 64 / 15.436677e9
target:    188 * 2.430e9 * 64 / 26e9
```

The checker reports both values.  A map-only construction above the target
allowance is not admitted to implementation by this frozen protocol, even
when every field operation, state access, branch, report, and loop instruction
is free.  This is a conditional screening rule, not a hardware lower bound.

The reference row demonstrates the model's transfer limit.  Its 2,276.25
charged source operations/update predict 12.844705 B/s at 29.23776 T
lane-operations/s, while the fused kernel measured 15.436677 B/s.  The
1.201793 ratio reflects compiler lowering and overlap that this static model
does not represent.  Exact field checks and source counts are proved below;
their conversion to throughput is explicitly modelled.

For the complete comparable B16 source ledger, the reference charges 1,230
word operations/update to selector plus `L_j` composition and 1,046.25 to the
current direct reductions: **2,276.25** in those stages.  The sparse candidate
charges 476.625 reducer operations and, with `PACKED_ALU_SQUARE=1`, 75 ALU
spread operations in addition to its synthesized maps.  Unchanged products,
inverse transforms, state, control, and reporting are omitted equally from
both rows; these are source-circuit counts, not SASS or throughput.

## Frozen matrices and oracle

The native C++17 checker reconstructs the field isomorphism from beta root
`0x30a16693fefe59e60962fc4e3ddc388eb`, then derives:

- sparse-to-normal, needed once/update for the invariant selector;
- sparse-to-beta and beta-to-sparse, each needed once per batch 16; and
- the eight direct sparse maps `L_j = I + sigma^j`, `j=3..10`, needed on both
  coordinates.

Every matrix is checked on all 131 input basis vectors and 4,096 deterministic
dense vectors against the repository arithmetic.  Rank and canonical output
bits are mandatory.

## Candidate circuit families

The following finite set is fixed before the canonical run:

1. **Generated masked diagonals.**  A nonzero `(diagonal, output-word)` term
   uses one shifted source word and one constant mask.  Within the checker's
   one-shifted-word abstraction, the direct count is one word shift plus one
   fused `LOP3(out, shifted, mask)` per term, except the zero diagonal needs no
   shift.  A 131-bit cross-word shift can cost more than that abstraction.
   The stronger oracle makes every shift free and charges one `LOP3` per term.
2. **Warp-divergent fixed circuits.**  With near-uniform `j`, the checker
   reports the expected number of distinct branches in a 32-lane warp and the
   serialized sum of the eight fixed circuits.  It also reports an impossible
   oracle control in which every update receives the cheapest `L_j` with free
   selection and compaction.
3. **All-output common-shift circuit.**  Shifted words are shared across all
   eight maps, exact equal masks are reused, all eight outputs are materialized,
   and a 5-word, eight-way branchless selector is charged.  The checker counts
   union terms, mask signatures, output memberships, and selector operations.
4. **Greedy XOR common-subexpression circuits.**  A deterministic Paar-style
   heuristic synthesizes sparse-to-normal, the cheapest fixed `L_j`, and their
   joint outputs.  Each bit-level XOR is optimistically charged as one full GPU
   ALU instruction and all intermediate storage/bit extraction is free.  This
   is a candidate construction, not a universal XOR lower bound.
5. **Frobenius composition identities.**  The checker verifies
   `L_(a+b)=L_a + sigma^a L_b`, `L_(2a)=L_a composed with L_a`, and
   `L_(j+1)=sigma L_j + L_1`, then prices sequential powers, a branchless
   `3 + 1 + 2 + 4` chain, and an impossible perfectly compacted mean-`j`
   control.  A sparse ALU square is charged as five 15-operation spreads plus
   the 82-operation two-fold reducer.

The checker also reports identical rows, identical diagonal masks, pairwise
matrix distances, and the live intermediate count.  This prevents a claimed
common-subexpression win from being inferred from matrix density alone.

## Decision and stop rule

A candidate is admitted to a separate GPU protocol only if all exact checks
pass and at least one concrete circuit satisfies both frozen model rules:

1. its complete selector/direct-map/conversion map cost is below the 26 B/s
   source-operation screening allowance; and
2. after adding the frozen 476.625 sparse-reducer and 75 ALU-square source
   operations, its comparable B16 subledger is below the reference's
   2,276.25 operations/update.

The first rule gives all other walk work zero cost but remains a modelled
admission rule rather than a universal necessary condition: actual compiler
lowering and pipeline overlap can differ by operation mix.  Passing authorizes
a separate implementation protocol, not a speed claim.  If every concrete
family fails, the workstream stops for the enumerated constructions with a
compact negative result and no GPU/full-walk experiment.
