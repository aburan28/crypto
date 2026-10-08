# P-256 variable-column S17 multilevel selector, round 19: result

Date run: 2026-10-05

**EXACTNESS PASS; ALL ATTACK-PROMOTION GATES FAIL.**  A global variable-column
selector is a real exponent improvement over round 18's fixed-atom scan, but
it is not close to an end-to-end P-256 attack.  The optimistic 138,031-row
projection falls from round 18's `2^272.404` field multiplications to
`2^165.156` P-256 field-multiplication equivalents (FME).  That is a
107.248-bit stage improvement, while remaining 45.156 bits above the frozen
`2^120` collection gate and 32.690 bits above the same-unit rho reference.

The proposed Dickson-aware 4+4+4+5 schedule does not beat direct 8+9.  Every
active leaf passes its depth-18 Dickson chain, so the chain removes **zero**
global combinations.  At `B=20` the structured arm performs 1.2391 times the
baseline's counted FME, writes 7.799 GiB, reads 9.030 GiB, and materialises a
4.515-GiB one-target frontier.  External bucketing lowers resident bucket
memory, but moves the same list to disk and adds work.  It is a relabelling,
not a structural advance.

## One boundary table, one unit

One FME is one P-256 field multiplication; a complete projective addition is
charged as 17 FME.  Sparse scalar operations receive an optimistic one-FME
minimum.  Sorting and byte I/O remain separate, so every projected FME total
below is a lower bound.  Rho is `1.3*sqrt(n)` group additions, converted to
the same unit.

| variant | status | projected total FME, log2 | log2(total / rho) | projected materialised bytes, log2 | correctness | class |
|:--|:--|--:|--:|--:|:--:|:--|
| Pollard rho reference | reference | 132.466 | 0 | small | generic reference | reference |
| round-18 fixed 16+1 atom scan | rejected | 272.404 | 139.938 | not dominant | bounded scan exact | accounting |
| exact balanced 8+9 | rejected | **165.156** | **32.690** | 134.092 left table | yes | stage advance versus fixed atoms |
| Dickson-aware 4+4+4+5 external join | rejected | **at least 165.156** | **at least 32.690** | **148.926 frontier** | yes | relabelling |

The common optimistic per-target lower bound is `2^148.028` FME.  Dividing by
the 0.9640001820745058 modeled relation probability gives `2^148.081` FME per
usable relation, 45.081 bits above the `2^103` gate.  The structured storage
projection is 98.926 bits above the `2^50` gate.  It fails even if right-sum
construction is shared across every collection target and byte I/O costs
nothing.

## Measured exact prefixes

All rows use actual P-256 points selected from `FB1h2f8621cda105`.  The
baseline holds the complete eight-sum index in memory.  The structured arm
builds exact signed 4- and 5-lists, joins 4+4 and 4+5, and processes all 256
x-coordinate buckets for each target.

| active B | signed 8-list | signed 9-list | baseline FME | structured FME | structured / baseline | baseline logical left bytes | structured materialised peak bytes | exact |
|--:|--:|--:|--:|--:|--:|--:|--:|:--:|
| 17 | 6,223,360 | 12,446,720 | 903,818,583 | 1,113,769,569 | 1.2323 | 273,827,840 | 765,473,280 | yes |
| 18 | 11,202,048 | 24,893,440 | 1,779,616,529 | 2,197,217,483 | 1.2347 | 492,890,112 | 1,479,915,008 | yes |
| 19 | 19,348,992 | 47,297,536 | 3,337,713,009 | 4,128,581,739 | 1.2369 | 851,355,648 | 2,732,507,648 | yes |
| 20 | 32,248,320 | 85,995,520 | 6,002,572,382 | 7,437,716,615 | **1.2391** | 1,418,926,080 | 4,847,997,440 | yes |

The measured `log2(work)` slopes against `log2(B)` are 11.6502 for direct
8+9 and 11.6841 for 4+4+4+5 over `B=17..20`.  The narrow-prefix fits are
diagnostics.  The challenge-size projection instead uses the exact widths

```text
C(131458,8) * 2^8 =
566139037710260120415636851490531741696       (2^128.734)

C(131458,9) * 2^9 =
16537550334891931739696769806317866099097600  (2^143.569).
```

The result JSON's `peak_resident_logical_bytes` field counts the simultaneously
loaded left/right bucket records.  The structured implementation also retains
its reusable projective 4- and 5-lists.  From their frozen record counts and
the source's 112-byte `Partial` layout, the B=20 total logical candidate peak
is 84,720,940 bytes (80.80 MiB), before allocator and factor-base overhead;
the materialised disk frontier remains 4,847,997,440 bytes.  This accounting
clarification does not affect either projected gate.

## Exactness and target outcomes

- Every baseline and structured cell returns the same relation set.
- `B=17` independently exhausts 131,072 complete signed 17-sums; `B=18`
  independently exhausts 2,359,296.  Both references equal both selectors.
- The planted witness is recovered and freshly replayed in all eight selector
  cells and both independent-reference cells.
- Every hash-public cell has zero relations.  That is a bounded observation,
  not a full-factor-base nonexistence claim.
- False positives / false negatives: `0 / 0` throughout.
- The planted relation digest is
  `4679da63b1051820593722c2d2fa0d99370aa60e9cd7ca990dee361af2ab20cc`;
  the empty public set has the standard SHA-256 empty digest.

The B=19 and B=20 cells do not have the separate 17-subset reference.  Their
completeness evidence is agreement between independently implemented direct
8+9 and multilevel 4+4+4+5 constructions, plus the exhaustive smaller cells,
exact combinatorial list-width checks, and fresh replay of every emitted hit.

## What the Dickson structure did

For every active factor-base x-coordinate the native harness recomputes

```text
D_2(x) = x^2 - 2
```

through 18 levels and checks the frozen terminal.  Column masks propagate
those certificates through every 4-, 5-, 8-, and 9-list record.  Survival is
1.0 at all four sizes: factor-base membership is a per-leaf constraint and
does not correlate the choice of 17 different columns.  Exact intermediate
P-256 points enforce the split summation-polynomial equations, but computing
them still emits every generic 8- and 9-sum.

The degree evidence therefore stays exactly where it was:

- local split join degree: 2;
- frozen structured residual maxima at depths 1, 2, 3: 3, 3, 4;
- unsplit S17 degree of regularity: unknown.

Degree at most 5 passes.  It does not compensate for an unfiltered
`2^143.569` right frontier.

## Relation matrix and decision

The 138,031-row sparse matrix still has 2,346,527 nonzeros and a 24.278-MiB
CSR-plus-three-vectors model.  The conservative Wiedemann schedule contributes
616,939,492,732 nonzero additions and Berlekamp--Massey contributes
17,281,205,764 scalar operations.  Even charging each at one FME leaves the
total at `2^165.156`; relation acquisition dominates.

Do not attempt the full-depth unplanted P-256 relation.  The exact variable-
column construction clears the fixed-atom obstruction but fails the
per-relation, collection, storage, and rho gates.  A credible successor needs
a non-separable cross-column invariant that rejects candidates before the
8/9 frontier is emitted.  Another bucket layout or another local Gröbner
degree reduction cannot change this verdict.

## Reproduction and host

Source commit used for the run:
`7536268c808cbfa7714ffca090a17b3e8315a6d2`.

```bash
cargo test --bin p256_s17_multilevel_selector
cargo clippy --bin p256_s17_multilevel_selector -- -D warnings \
  -A clippy::mismatched_bit_width_type
cargo build --release --bin p256_s17_multilevel_selector
target/release/p256_s17_multilevel_selector \
  --round18 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_outer_scan_round18_20261005/outer-result.json \
  --round6 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_residual_scaling_round6_20261004/degree-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_multilevel_selector_round19_20261005/selector-result.json \
  --temp-parent /tmp
```

The count run used Rust 1.99.0 (`b940084d7`, LLVM 23.1.1), Linux x86-64,
five exposed cores of an AMD EPYC 9V74 VM, and 17 GiB RAM.  No wall-time claim
is made.  `selector-result.json` is 22,056 bytes with SHA-256
`3096540621408e4a48cfa18963ad01da9686d3527ee26776c8b6cf6f45a71114`.
