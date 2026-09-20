# The Riemann–Roch encodings in front of WDSat

**Exporter:** `ecc2k130/codegen/wdsat.py`, tests `ecc2k130/codegen/testwdsat.py`
**Experiment:** `ecc2k130/codegen/wdsat_rr_panel.py`
**Frozen evidence:** `research/wdsat_rr_20260920/`
**Background:** [`research/nagao_relations/README.md`](../../nagao_relations/README.md)
(the incidence and norm encodings, and the panel that did not complete),
[`RESEARCH_ECC2K130_RR_SOLVER_PANEL.md`](RESEARCH_ECC2K130_RR_SOLVER_PANEL.md)
(the other Riemann–Roch line — a search procedure, not an encoding),
`research/index_calculus_baseline_20260914/regression/README.md`
(the frozen WDSat suite this round is *not* a run of)

**The question.** WDSat is a SAT solver written for Weil-descent instances;
this repository vendors and builds it, but only ever fed it the S4 corpus.
The Riemann–Roch frontends are Boolean systems of exactly the shape it was
built for, and one of them is recorded as unfinished: in the 2026-09-13
comparison **all 24 nine-bit RR-incidence instances hit the ten-second
CryptoMiniSat watchdog**, while RR-norm completed 110/110. That is an open
question about the solver, not about the encoding, and it is the one this
round can settle.

## 1. Boundary and unit

*Unit.* **Conflicts**, each solver's own counter — the same unit the frozen
S4 regression compares, and the only one both engines report natively.
Wall-clock seconds are recorded alongside as a practicality note
(`AGENTS.md` §6), never as the metric.

*Floor.* Exhaustive search over the free Boolean variables: `2^v`
assignments for a system with `v` free bits. No solver deciding this
system can be credited below it, and it is what makes "the solver is doing
algebra rather than search" checkable.

*Reference.* **CryptoMiniSat on the identical equations**, through the same
`ir.Prog` the WDSat arm exports, with the same target, formulation and
budget. This is the comparator the RR work already uses.

## 2. What could not be matched, and why

WDSat's ANF mode — `__XG_ENHANCED__`, which is what makes it worth using
and what the frozen build defines — **parses no OR-clauses**. In
`dimacs.c`, a line beginning `x` is an XOR of monomials and *every other
line is read as a single-literal unit clause*; there is no branch for a
general disjunction.

The RR encoding's `encode()` adds four things that are disjunctions: the
Hamming-weight bound `atMost(v, weight)`, `b ≠ 0`, the distinct-and-nonzero
abscissa domain, and the lexicographic ordering that breaks the summand
symmetry. **None of them can be given to WDSat in this mode.**

So the matched panel compares the two solvers on the **equational part
only**, and the models are filtered afterwards. The domain constraints are
a real advantage that CryptoMiniSat has and WDSat structurally cannot be
given here, and §5 reports the CNF arm with them as a separate,
differently-labelled row rather than folding it into the comparison.

This is also why this round is **not** a run of the frozen S4 regression.
That suite pins its inputs, its binary and its ANF semantics; a different
encoding is a different instance family and, because WDSat's structures
are statically allocated, a differently compiled binary. `AGENTS.md` §8
anticipates exactly this and asks for an equivalent matched suite instead,
which `research/wdsat_rr_20260920/` is.

## 3. Falsification target, declared before measuring

The round is a **success** if, at the ten-second budget the earlier
comparison used:

* WDSat **decides** — SAT with a verified witness, or UNSAT — at least one
  RR-incidence instance at `n = 9` that CryptoMiniSat leaves undecided on
  the identical equations; and
* every model WDSat returns reconstructs a genuine relation through
  `nagaodecomp.reconstructWitness`, with zero invalid witnesses admitted.

It is **abandoned** if WDSat leaves undecided every instance
CryptoMiniSat leaves undecided, under both the Gaussian-elimination-on and
-off configurations the frozen suite uses.

Inadmissible: changing the encoding, target, formulation or budget between
arms; counting a timeout as a refutation; reporting the CNF arm's extra
domain constraints as part of the solver comparison; and any claim about
`2^131`, since this prices a decomposition oracle and relation collection
stays `Θ(2^n)` with any oracle polynomial in the factor base.

## 4. Results

`python3 wdsat_rr_panel.py --n 9 --targets 8 --budget 10`, eight targets per
formulation, both arms on the identical equations.

| n | formulation | arm | decided / 8 | witnesses verified | bad | conflicts, median |
|--:|:--|:--|--:|--:|--:|--:|
| 9 | norm | WDSat, `-x` (Gaussian) | **8** | 8 | 0 | **62** |
| 9 | norm | WDSat, plain | **8** | 8 | 0 | 79,151 |
| 9 | norm | CryptoMiniSat | 8 | 8 | 0 | not reported |
| 9 | incidence | WDSat, `-x` | **0** | — | — | did not run (§4.2) |
| 9 | incidence | WDSat, plain | **0** | — | — | did not run (§4.2) |
| 9 | incidence | CryptoMiniSat | 7 | 7 | 0 | not reported |

At `n = 5` every arm decides every instance and every witness verifies,
4/4 on both formulations — the sanity rung, in the frozen output.

### 4.1 Gaussian elimination is worth three orders of magnitude

The one clean positive. On RR-norm at `n = 9`, XORGAUSS cuts the median
conflict count from **79,151 to 62**, a factor of **1,277**; at `n = 5` the
same comparison is `488 → 18`, a factor of 27. Both arms are the same
binary on the same instance with `-x` on and off, which is exactly the
paired comparison the frozen S4 suite makes, so the measurement means
here what it means there. This is the property WDSat exists for, and the
Riemann–Roch norm encoding has it.

It does not make WDSat the faster engine end to end: CryptoMiniSat decides
the same instances in a median 0.33 s against WDSat's 2.04 s. Conflicts
and seconds disagree because they measure different things, which is
§4.3.

### 4.2 RR-incidence does not fit the solver at all

Not a timeout and not a refutation: the process dies with SIGSEGV before
printing anything. The cause is WDSat's static allocation, and it is
measurable independently of any encoding. `monomials_to_column` is
`__MAX_ANF_ID__ × (__MAX_ID__ + 1) × (__MAX_DEGREE__ − 1)` of
`uint_fast64_t`, and a one-variable trivially satisfiable probe built at a
range of sizes (`ceiling.txt`) puts the boundary between **22.9 and 23.9
million entries**, about 370 MB:

| `__MAX_ANF_ID__` | `__MAX_ID__` | product (M) | result |
|---:|---:|---:|:--|
| 3,938 | 5,314 | 20.9 | ok |
| 4,300 | 5,314 | 22.9 | ok |
| 4,500 | 5,314 | 23.9 | SIGSEGV |
| 4,191 | 5,729 | 24.0 | SIGSEGV |
| 4,875 | 6,656 | 32.5 | SIGSEGV |

RR-norm at `n = 9` needs 3,937 ANF variables and 20.9 M entries, and fits.
RR-incidence needs 4,874 and 32.5 M, and does not. Dropping the three
`p_i ≠ x_R` witnesses (`--core-witnesses`, filtered from models
afterwards) gets incidence to 4,190 variables — under the variable
ceiling, still over the product ceiling at 24.0 M, still SIGSEGV.

So **the pre-registered success condition cannot be met for incidence**,
and not because incidence is hard: it does not reach the solver. The
condition asked for WDSat to decide an incidence instance CryptoMiniSat
leaves undecided. In this panel CryptoMiniSat decides 7 of 8, and WDSat
decides none. The round is **abandoned on its own terms** for incidence
and reports a structural boundary instead.

That boundary is the reusable result. Any encoding wanting WDSat must fit
about 23 M `monomials_to_column` entries, which for a quadratic ANF is
roughly `4,300` unary variables — a hard budget that the Weil-restricted
ONB multiply spends quickly.

### 4.3 The declared unit is not available on the reference arm

§1 named conflicts, because it is the frozen suite's unit and both engines
were expected to report it. **They do not.** WDSat prints its counter on
every run; CryptoMiniSat through `pysat`'s `accum_stats` returns no
conflict figure for this solver, so the column reads "not reported" above
and the cross-solver comparison **cannot be made in the declared unit at
all**.

What survives is the within-WDSat comparison of §4.1, which is a paired
measurement in the declared unit and is the only conflict number this note
claims. Seconds are reported for both arms as a practicality note and are
not promoted to the metric to cover the gap; doing that would be changing
the unit after seeing the result, which `AGENTS.md` §6 exists to prevent.

### 4.4 Two exporter bugs the witness gate caught

Both would have produced a solver that ran, answered, and meant something
else:

* **The parity convention is inverted.** WDSat's ANF equations assert that
  the listed terms XOR to *one*; `dimacs.c` absorbs a `T` by negating the
  first literal. An exporter written for "XOR to zero" produces a
  satisfiable instance describing the complementary system.
* **`ir.Prog.addInput` does not deduplicate.** One input reference can
  name several nodes, and `emitCnf` collapses them through its literal
  map. The first exporter gave the second node its own variable, so the
  domain conditions written about `p_0` constrained a different `p_0`.
  Every equation was satisfied and every witness was junk.

Neither is visible from the solver's output. They are why `testwdsat.py`
feeds models back through the original `ir.Prog` rather than trusting a
SAT answer, and why the panel re-derives each relation with
`nagaodecomp.reconstructWitness`.

## 5. The CNF arm with its domain constraints

Reported separately, because it is not part of the solver comparison.
CryptoMiniSat can additionally take the disjunctive domain the CNF
frontend builds — the Hamming-weight bound, `b ≠ 0`, distinct abscissas
and the lexicographic symmetry break — none of which WDSat's ANF mode can
read. Withholding them from both arms is what makes §4 matched; it also
makes both arms solve a weaker system than the RR frontend normally does.

With the equational domain of §2 in place, the filter still rejects
models: at `n = 9` under `--core-witnesses`, CryptoMiniSat returns 8 of 8
decided on incidence but only 2 of 8 witnesses reconstruct, the other 6
being models with `p_i = x_R` that the ANF no longer forbids. Those are
counted as bad witnesses rather than quietly dropped. The full witness set
restores 7 of 7.

## 6. Classification

By `AGENTS.md` §3:

| result | what moved | class |
|:--|:--|:--|
| XORGAUSS worth 1,277× conflicts on RR-norm | a solver constant, within one engine | **engineering** |
| RR-incidence unreachable at `n = 9` | nothing; a tool boundary was measured | **accounting** |
| The declared unit unavailable on the reference arm | the comparison was not made | **accounting** — stated, not worked around |
| Anything about `2^131` | nothing | not measured; §1 |

No row is an advance. Relation collection stays `Θ(2^n)` with any oracle
polynomial in the factor base, and this prices a decomposition oracle.

## 7. What would make this worth another round

* **Shrink the encoding under 23 M entries.** The ONB multiply dominates
  the variable count; a subspace factor base of dimension `l` would
  replace `m`-bit abscissas with `l`-bit ones, which is what the frozen S4
  corpus does and why it fits.
* **A conflict counter for the reference arm**, so §4.3's comparison can
  actually be made. `pysat` is the wrong binding for it; the frozen
  suite's own runner reads WDSat's counter directly and would need a
  CryptoMiniSat equivalent.
* **`__FIND_ALL_SOLUTIONS__`**, which the vendored `config.h` leaves off.
  A decomposition oracle wants every root, not the first, and the panel
  currently stops at one.
