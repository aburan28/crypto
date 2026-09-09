**Full-orbit factor-base search — degree 19**

This continuation expands the four-orbit search from the earlier fixed pool of 37 signed Frobenius orbits to all 6,909 nonzero signed Frobenius orbits in the prime-order subgroup of

\[
K_1:y^2+xy=x^3+x^2+1
\quad\text{over}\quad
\mathbb F_2[z]/(z^{19}+z^5+z^2+z+1).
\]

It finds one 152-point quotient base that improves both the exactly-three relation support used by a pure symmetric `S_4` collector and the at-most-three support used when a separate two-summand branch is retained. These are exact finite results on public subgroup labels. No unknown scalar is recovered.

**Selected base.** The four public scalar-orbit representatives are `[203, 2143, 2901, 9853]`. Canonical point representatives, with integers written as polynomial-basis masks, are:

| scalar representative | u-coordinate | y-coordinate |
|---:|---:|---:|
| 203 | 16795 | 144921 |
| 2143 | 1315 | 90425 |
| 2901 | 8461 | 38471 |
| 9853 | 6685 | 369649 |

Close these points under Frobenius and negation. The resulting base has 152 nonzero points, eight ordinary Frobenius orbits, and four signed Frobenius relation columns. Its two-torsion-closed affine lift has 305 points and the same four projected columns. The implementation exposes an explicit-orbit constructor, records `point = (-1)^s*pi^k(representative)` for each member, uses `(-1)^s*lambda^k` in relation rows, and accepts caller-supplied nonlinear bases in the index-calculus driver.

| Four-orbit construction | Exactly three target orbits | At most three target orbits |
|---|---:|---:|
| Earlier unique optimum in the fixed 37-orbit pool | 6224 / 6909 (90.0853959%) | 6257 / 6909 (90.5630337%) |
| Full-universe selected base | **6256 / 6909 (90.5485598%)** | **6300 / 6909 (91.1854103%)** |
| Improvement at unchanged base size and column count | +32 orbits / +1,216 nonzero targets | +43 orbits / +1,634 nonzero targets |

Exactly 44 target orbits are covered by the one/two-summand branches but not by exactly three summands. The group-law oracle independently obtains the same two support counts and checks the complete symmetric-cube representation total.

**Search method and guarantee.** On the order-262543 subgroup, Frobenius acts as multiplication by the public characteristic-polynomial root `lambda`. Each candidate orbit therefore contains the 38 public labels `+-lambda^j*k`. For an orbit type `(i,j,k)`, normalize the first summand to its public representative and enumerate the other two 38-element orbits. Mapping each nonzero modular sum back to its signed Frobenius orbit gives the exact support of that type. A four-orbit base is the union of its 20 triple-type supports, plus its single and pair types for the hybrid objective.

For every current base and each of its four removal positions, the implementation scores all 6,906 candidates not already among the three retained orbits. It accepts the steepest improving exchange and repeats. Each terminal pass therefore covers all 27,620 proper one-orbit neighbors, plus four reinsertions of the removed current orbit. The terminal scans show that no one-orbit replacement improves the selected base under either objective. This proves one-exchange local optimality over the full candidate universe; it is not a global optimum over all `binomial(6909,4)` bases.

The search used deterministic two-orbit perturbations of the fixed-pool incumbent. The retained census contains nine exactly-three runs and seventeen hybrid runs. Seed 9118 is the unique record in that census for both objectives, and both objective-specific climbs terminate at the same base. The restart census is in `restart-summary.json`; complete JSONL outputs preserve every accepted exchange and terminal full-neighborhood scan.

**Membership polynomial.** Let `M_i(U)` be the degree-19 minimal polynomial of each selected u-coordinate and let `g(U)` be their product.

| u-coordinate | Minimal-polynomial mask | Descended Boolean-degree upper bound |
|---:|---|---:|
| 16795 | `0xd0aa9` | 3 |
| 1315 | `0xe7049` | 3 |
| 8461 | `0xf7761` | 3 |
| 6685 | `0xcf9b9` | 4 |

The squarefree product has degree 76 and mask

```text
0x11c1371440e1f138d579
```

It has exactly the 76 selected abscissae as roots. Its descended syntactic Boolean-degree bound is five, compared with six for the earlier fixed-pool base. This is a property of the membership equation, not a solving-degree estimate or a measured SAT speedup. The independent field-arithmetic certificate also verifies both rational-fiber trace predicates for every root.

**Independent checks.** `verify_full_orbit_result.py` reconstructs all 6,909 target orbits from the public characteristic polynomial, recomputes both supports using an independent modular enumeration, compares the two bitsets bit for bit, checks the group-law totals, and reads back all four terminal removal scans for each objective. `full_orbit_predicate_certificate.py` independently verifies the curve points, minimal polynomials, complete root set, trace predicates, squarefreeness, and Boolean-degree bound.

The selected explicit-orbit constructor also passes a deliberate native-XOR `S_4` round trip at degree 19. For the planted public target `[203+2143+2901]G`, the exact-domain SAT path returned a valid three-point decomposition in one solver call and one model, with 5,472 conflicts, no invalid model, and no exhausted budget. The test body took 2.91 seconds in the recorded run. This establishes operational coupling of the selected domain to the SAT encoder on one controlled instance; it is not a comparative SAT benchmark or a general runtime claim.

Reproduction commands:

```sh
cp research/sat_factor_base_review_20260908/continuation-01/source_snapshots/Cargo.lock Cargo.lock
cargo build --release --example koblitz_affine_quotient_base --locked
target/release/examples/koblitz_affine_quotient_base 19 --full-orbit-search research/sat_factor_base_review_20260908/continuation-01 6 --kick-seed 9118
target/release/examples/koblitz_affine_quotient_base 19 --full-orbit-search research/sat_factor_base_review_20260908/continuation-01 6 --exactly-three --kick-seed 9118
python3 research/sat_factor_base_review_20260908/continuation-02/full_orbit_predicate_certificate.py research/sat_factor_base_review_20260908/continuation-02/full-orbit-search-hybrid-kick-9118.jsonl
python3 research/sat_factor_base_review_20260908/continuation-02/verify_full_orbit_result.py
python3 research/sat_factor_base_review_20260908/continuation-02/summarize_restarts.py
cargo test --release cryptanalysis::koblitz_index_calculus::tests::selected_explicit_orbit_base_sat_round_trip_n19 --locked -- --ignored --test-threads=1 --nocapture
```

The final source-bound repetitions took 8.34 seconds for the hybrid objective and 7.83 seconds for the exactly-three objective on the recorded host. These timings describe the finite factor-base search implementation. They are not timings for solving the coupled Semaev SAT system.
