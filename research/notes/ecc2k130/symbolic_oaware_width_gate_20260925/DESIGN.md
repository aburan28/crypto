# Smallest complete symbolic relation and current interface gaps

## Full-point relation to encode

Over `F=GF(2^n)`, write an intermediate point as `(o,x,y)`. `o=1`
is canonical O with `x=y=0`; `o=0` must satisfy
`y²+xy=x³+a*x²+b`. Each factor is affine, with its x coordinate the
rotated slot's Boolean-linear basis combination; keep both rational y
signs through an unconstrained field y subject to the curve equation.
Set the initial prefix to O. For every addition `R=P+Q`, Q is an affine
factor; one of the following disjoint cases must hold:

| Case | Guard | Required output |
|:--|:--|:--|
| Identity | `oP=1` | `R=Q` |
| Inverse | `oP=0`, `d=xP+xQ=0`, `s=yP+yQ=xP` | `R=O` |
| Doubling | `oP=0`, `d=0`, `s=0`, `xP!=0` | unique `λ` with `xP*λ=xP²+yP`; `xR=λ²+λ+a`, `yR=xP²+(λ+1)*xR` |
| Distinct x | `oP=0`, `d!=0` | unique `λ` with `d*λ=s`; `xR=λ²+λ+d+a`, `yR=λ*(xP+xR)+xR+yP` |

For same x, the curve equation implies `s=0` or `s=xP`; when `xP=0`
these collapse to the inverse-to-O case, so division by zero never enters
doubling. A circuit must encode the branch selector, zero/nonzero guards,
field equations, canonical O code and all output equations bidirectionally.
The existence of these mathematical relations is **not** an implemented
complete circuit or a proof that a chosen Tseitin/CNF encoding is sound.
The terminal `(o,x,y)` equals the exact full Q+T target, not just its x.
Independent group-law replay must lift every SAT model and reconstruct each
of the 32 #781/#785 n13-m5 statuses. Negative statuses require externally
checked certificates; projected Q negativity needs all four T branches.

This affine route uses no explicit prefix-state enumeration. It does use
symbolic prefix variables. A chain with m factors has at least `m-2`
intermediate x fields before y/O/slope/selector and solver auxiliaries.
At n131, unequal m9 has `5*15+4*14=131` factor bits plus `7*131=917`
prefix x bits, a 1,048-bit raw-chain floor. Unequal m10 has
`14+9*13=131` factor bits plus `8*131=1,048` prefix x bits, a 1,179-bit
floor. Balanced m10 has 130+1,048=1,178. These floors concern this affine
chain layout; a direct non-chain exporter has a different state layout.

## What the current code actually supports

| Interface | Source | Current limit / semantic gap |
|:--|:--|:--|
| Binary curve arithmetic | `src/binary_ecc/f2m.rs`, `curve.rs` | multiword n131 and full O-aware point arithmetic exist numerically, not symbolic |
| Fast GF2 used by Semaev | `semaev_decomp.rs` | element/modulus are `u64`; `Gf2::new` asserts degree <=63 |
| Symbolic field polynomials | `pq_descent_symbolic.rs` | coefficient and monomial masks are `u64`; max 64 Boolean vars, only 2/3 summands; S3/S4 x-only |
| Boolean Gröbner monomials | `pq_groebner_f2.rs` | `F2BoolMono.mask:u64`, variable index `<64` |
| Boolean-system SAT wrapper | `semaev_sat.rs` | `n_vars<=64`, monomial key and model assignment `u64` |
| CDCL core | `sat.rs` | `u32` variable IDs/i32 literals can address more than 64; wrapper/exporter is the blocker |
| #781 O-aware n13 exporter | `rotated_s3_o_branch_20260925/export.py` | branch-complete on toy but explicitly enumerates rational factor x and reachable prefix x/O states |
| #785 32-target panel | `n13_oaware_sat_benchmark_20260925` | fixed complete toy statuses/model checks; no symbolic n131 exporter |

The narrower first engineering increment is a **solver-neutral XOR/AND
expression DAG** with `u32` node/variable IDs, incremental CNF or native-XOR
emission, and multiword GF(2^n) bit-vector operations (linear
squaring/reduction and quadratic multiplication). It can bypass the old
`F2BoolMono`/F4 polynomial representation rather than widening that whole
engine at once. A separate, larger route would widen `Gf2`, every Boolean
monomial container, model assignment and blocked-model API before reusing
the existing ANF/Groebner/SAT wrapper; porting only a `u64` mask type is
insufficient. The DAG route still needs the full-point guarded addition
relation above, DIMACS/Tseitin emission, and a **fail-closed model-to-point
map**: reject missing/unassigned factor bits, noncanonical O triples,
off-curve y values, unsatisfied branch guards, inconsistent prefix points,
and terminal point/sign mismatch; then replay all factors with independent
group arithmetic. Every
guard needs truth-table tests on tiny complete curves before the #781/#785
same-input 32-target comparison. Admission then needs export size/time/RSS
caps and independent model and negative-proof replay; no estimate of n131
SAT success follows from passing the width gate.

## μ4 alternative and its limit

`coordinate_search.rs` and `RESEARCH_EXOTIC_COORDINATES.md` already
investigate Kohel μ4 as a coordinate/tuple idea. The exact split μ4 map in
the protocol is useful because K0 has rational order-4 point `(1,0)`.
With `X2=1`, both quadrics reduce exactly to K0; `X2=0` is solely O.
This proves a full-point birational chart and the frozen toy exhaustive
check tests implementation details. Kohel's bidegree-(2,2) addition-law
*basis* has individual exceptional divisors; a circuit needs a complete
selector/coverage proof and projective scaling constraints. Binary Edwards
complete-formula claims concern a different model and likewise require
an exact K0 point conversion and Boolean-circuit validation. Neither is
currently the smallest **implemented** path. The affine relation above is
the smallest specified change that reuses K0's existing coordinate model.

Primary sources: [Kohel, *Efficient arithmetic on elliptic curves in characteristic 2*,
2016](https://arxiv.org/abs/1601.03669), especially Theorems 6.1 and 7.2;
[Bernstein et al., *Binary Edwards Curves*, 2008]
(https://cr.yp.to/newelliptic/edwards2-20080611.pdf). Existing repository
coordinate-search work predates this gate; no novelty is claimed for the
normal form itself.
