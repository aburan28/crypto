# Different-j and rational-fiber proof map

The native enumeration in `FIBER_AUDIT_V1.json` checks this exact map over
`F_(67^2)`. The Lean file in `lean/DifferentJ67.lean` uses the pinned Lean
4.34.1 toolchain and no external dependencies. These declarations prove
finite arithmetic and conditional scalar facts; they are not a formalized
elliptic curve or Vélu isogeny.

| Written statement | Lean declaration | Exact status |
| --- | --- | --- |
| Frozen source `j=7` differs from destination `j=11` | `DifferentJ67.distinctJ` | Checked residue inequality; missing formal bridge from the concrete curve formulas to the j residues and the general j-isomorphism theorem |
| Degree 3 and subgroup order 61 are coprime | `DifferentJ67.targetCoprimeKernel` | Fully checked arithmetic |
| A subgroup scalar killed by 3 is zero | `DifferentJ67.kernelScalarVanishes` | Checked in `Fin 61`; the identification of this scalar ring with the concrete curve subgroup is a written, unformalized bridge |
| Degree-three kernel collisions imply injectivity on the order-61 scalar subgroup | `DifferentJ67.injectiveOfThreeKernel` | Checked conditional structural implication: if `image a = image b` implies `3*(a-b)=0` in `Fin 61`, then the image function is injective. The kernel implication and concrete image function for the Vélu map remain unformalized. |
| Enumerated source count is three times the image count | `DifferentJ67.recordedFiberCardinality` | Checked arithmetic only; the 4575-point enumeration and 1525 fiber keys are independently executed native facts, not certified by this Lean theorem |
| Group order factors as `75*61` | `DifferentJ67.recordedSubgroupOrder` | Checked arithmetic only; group-order enumeration remains native evidence |

The central unformalized theorem would define the concrete `Fp2` field,
short-Weierstrass group laws, the rational three-point kernel and explicit
Vélu map, prove its homomorphism and kernel exactly, and then derive equal
cardinality of every nonempty fiber and injectivity on the order-61 subgroup.
The finite native audit checks the conclusion on this one field, including
all 3721 ordered subgroup pairs. It does not supply a general Lean theorem.

`lake build DifferentJ67` passed in `lean/BUILD_V4.stdout` under pinned Lean
4.34.1. `#print axioms` reports only standard `propext` and `Quot.sound`
for the new conditional theorem, with no `sorryAx` or custom claim axiom.
Earlier `BUILD_V2` and `BUILD_V3` outputs retain two failed proof attempts
that used identifiers unavailable under this pinned `Init` import.
`lean4checker` was unavailable on this host; the Lean kernel build is the
recorded check.
