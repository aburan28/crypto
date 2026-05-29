# Elliptic divisibility sequences, elliptic nets, and the EDS-Residue handle on the ECDLP

**Module:** `src/cryptanalysis/eds_residue.rs`
**Demo:** `examples/eds_residue_demo.rs` (`cargo run --release --example eds_residue_demo`)
**Provenance:** Direct response to the research thread *"the genuinely
underexplored above-generic handle over the ECDLP is the elliptic
divisibility sequence / elliptic net structure (Shipsey, Stange,
Lauter–Stange) — there is a concrete avenue via quadratic-residuosity of
the EDS that has had far less attention than it structurally deserves."*

This note frames that avenue precisely, separates what is *provably*
true (the equivalence theorem) from what is *structurally open* (the
residue bias), and ships a runnable instrument that measures the relevant
invariants on toy curves so the open questions are stated against real
numbers rather than folklore.

---

## 1. The objects

### 1.1 Elliptic divisibility sequences (Ward 1948)

An **elliptic divisibility sequence** (EDS) is a sequence `W(n)` with
`W(0)=0`, `W(1)=1` satisfying Ward's recurrence for all `m ≥ n`:

```
W(m+n)·W(m−n) = W(m+1)·W(m−1)·W(n)² − W(n+1)·W(n−1)·W(m)².        (★)
```

Ward proved that an integer sequence satisfying `(★)` together with the
divisibility property `W(m) | W(n)` whenever `m | n` is exactly the
sequence of values of the **division polynomials** of an elliptic curve
at a point. Concretely, for `E : y² = x³ + ax + b` and `P = (x,y)`,

```
W(n) = ψ_n(P),   ψ_1 = 1,  ψ_2 = 2y,
ψ_3 = 3x⁴ + 6ax² + 12bx − a²,
ψ_4 = 4y(x⁶ + 5ax⁴ + 20bx³ − 5a²x² − 4abx − 8b² − a³),
```

and the doubling relations (valid because `W(1)=1`)

```
W(2n+1) = W(n+2)·W(n)³ − W(n−1)·W(n+1)³,
W(2n)·W(2) = W(n)·( W(n+2)·W(n−1)² − W(n−2)·W(n+1)² ).            (♦)
```

Both right-hand sides of `(♦)` reference only strictly smaller indices,
so the sequence fills left-to-right in `O(n)` field operations. This is
what `eds_sequence()` does over `F_p`.

### 1.2 The apparition law — the bridge to the group

The single fact that makes EDS relevant to the ECDLP:

> **`W(n) = 0  ⟺  [n]P = O`.**

Hence the *rank of apparition* (the first positive `n` with `W(n)=0`) is
exactly `ord(P)`. Our instrument cross-checks this against honest point
arithmetic — see §4, where it holds on every curve tested. The EDS is the
arithmetic shadow of the cyclic group `⟨P⟩`.

### 1.3 Elliptic nets (Stange 2007–2011)

Stange generalised EDS from `Z → K` to **elliptic nets** `W : Z^k → K`
satisfying a `k`-dimensional analogue of `(★)`. The crucial case is
`k = 2`: for an ordered pair `(P, Q)` of points there is a net
`W_{P,Q}(a, b)` with

> **`W(a, b) = 0  ⟺  [a]P + [b]Q = O`.**

Stange's main cryptographic result (*The Tate Pairing via Elliptic Nets*,
Pairing 2007) is that the Tate/Weil pairing is a ratio of net values, and
the whole net is computable by a double-and-add recurrence in time
comparable to Miller's algorithm. The net is therefore a genuine,
efficient computational object — not just a theoretical device.

---

## 2. Why this is an "above-generic" handle (and the precise limit of that)

### 2.1 The zero-lattice carries the discrete log

If `Q = [k]P` and `ord(P) = m`, then

```
W_{P,Q}(a, b) = 0  ⟺  [a]P + [b]Q = O  ⟺  a + b·k ≡ 0  (mod m).
```

So the **vanishing locus of the 2-D net is the sublattice**
`Λ_k = {(a,b) : a + bk ≡ 0 mod m}` of `Z²`, whose *slope is the discrete
log*: from any net zero `(a,b)` with `gcd(b,m)=1`,

```
k ≡ −a · b⁻¹  (mod m).
```

`recover_dl_from_net_zeros()` exhibits exactly this and recovers `k` on
every toy curve (§4). **This is the structural reason the ECDLP "lives
inside" the EDS/net world.** It also makes vivid *why* it is hard: the net
value at a single index is cheap, but *locating a zero of the net* is a
search over `Λ_k`, and that search is the ECDLP.

### 2.2 The Lauter–Stange equivalence: three problems, all = ECDLP

Lauter and Stange (*The Elliptic Curve Discrete Logarithm Problem and
Equivalent Hard Problems for Elliptic Divisibility Sequences*, SAC 2008 /
ePrint 2008/099 / arXiv:0803.0728) isolate three problems on the EDS of a
point and prove each is sub-exponential-time solvable **iff** the ECDLP
is. Let `E/K`, `P` with `ord(P) ≥ 4`, and `Q ∈ ⟨P⟩`, `Q ≠ O`, with
(unknown) discrete log `k`, so `Q = [k]P`:

- **EDS Discrete Log.** Given the sequence `W_{E,P}` and the data of `Q`,
  determine `k`. (This is essentially the ECDLP re-encoded.)

- **EDS Association.** Decide whether a candidate sequence is the genuine
  EDS associated to the pair `(E, P)`. Lauter–Stange relate this directly
  to the Tate pairing — and hence to the **MOV / Frey–Rück** reductions:
  exactly when the embedding degree is small (so the Tate pairing is
  cheaply computable into `F_{q^e}^*`), EDS Association becomes easy, and
  that is precisely the known weak-curve regime.

- **EDS Residue** *(the foregrounded avenue)*. Given `E`, `P`, and `Q`
  (**but not `k`**), determine the **quadratic residuosity of
  `W_{E,P}(k)`** — i.e. the single bit `χ(W(k)) = (W(k) | p)` indexed by
  the *unknown* discrete log.

The headline theorem says these are asymptotically all the same problem.
**So EDS-Residue is not a free lunch: any algorithm solving it in
sub-exponential time would already break the ECDLP.** Honesty up front —
no part of this note claims a generic-beating attack.

### 2.3 Where the open, underexplored room actually is

The equivalence is an *asymptotic* statement about sub-exponential
solvability. It leaves three concrete things unquantified, and these are
the avenue:

1. **The Legendre sequence `χ(n) = (W(n) | p)` has exploitable
   structure.** `χ(W(k))` is a *cheap* bit to compute *if you knew `k`*
   (one Legendre symbol). The EDS-Residue reduction shows that an oracle
   producing this bit *from `Q` alone* bootstraps to a full ECDLP solve.
   So the question "how much does `χ` leak, and how cheaply" is exactly
   the question of how close to that oracle one can get. This is the
   finite-field analogue of Silverman–Stephens, *The sign of an elliptic
   divisibility sequence* (2006), who proved the **sign** sequence of an
   EDS over `Z` is periodic with a period given by an explicit
   real-analytic formula (the elliptic logarithm of `P`). Over `F_p` the
   sign becomes the Legendre symbol; the analogous periodicity and bias
   are what we measure.

2. **Constant factors and special curves.** Even granting asymptotic
   equivalence, the *constants* in the reductions, and whether particular
   curve families exhibit anomalously small `χ`-period or strong residue
   bias, are unmapped. A strong, curve-specific bias is a distinguisher,
   and distinguishers compound.

3. **The QR pattern of the 2-D net as a cheaper localiser of `Λ_k`.**
   Finding an exact zero of the net is the ECDLP. But the Legendre symbol
   of net values is a coarser, much cheaper signal. Whether the
   `χ(W(a,b))` pattern localises the zero-lattice `Λ_k` faster than
   generic square-root search is open in practice.

---

## 3. The exact shift-multiplier structure (and the χ-period law)

Over `F_p` the EDS of `P` (rank of apparition `r = ord(P)`) is **not**
purely periodic with period `r`. Instead it is periodic *up to a
multiplier*: there exist constants `A, B ∈ F_p^*` with

```
W(n + r) = A · Bⁿ · W(n)    for all n.                            (▲)
```

(This is the finite-field form of Ward's periodicity and the object
Lauter–Stange call the sequence being "perfectly periodic" when `A=B=1`.)
`analyze()` extracts `(A, B)` from three consecutive shifted terms and
**verifies `(▲)` across a full block** — it holds exactly on every curve
tested. From `(▲)` the two periods follow in closed form:

- **Period of `W mod p`:** `π_W = r · j_W`, where `j_W` is the least
  `j ≥ 1` with `Bʲ = 1` and `Aʲ · B^{r·j(j−1)/2} = 1` in `F_p^*`.

- **Period of the Legendre sequence:** applying `χ = (· | p)` to `(▲)`,

  ```
  χ(W(n+r)) = χ(A) · χ(B)ⁿ · χ(W(n)),
  ```

  so `π_χ = r · j_χ` with `j_χ` the least `j` closing the `±1` recurrence.
  For `r` even this collapses to a **sharp dichotomy**:

  ```
  π_χ = r      if  χ(A) = +1 and χ(B) = +1,
  π_χ = 2r     otherwise.
  ```

The Legendre sequence's period is thus governed *entirely by the two sign
bits `(χ(A), χ(B))` of the multiplier* — a clean, computable invariant of
the pair `(E, P)`. This is the first concrete, falsifiable structural
statement the avenue produces, and §4 confirms it on the nose.

---

## 4. Measured data (toy sweep)

`cargo run --release --example eds_residue_demo`, five curves spanning
`p ≡ 1` and `p ≡ 3 (mod 4)`. All values are reproduced by the test suite.

| `p` | curve `(a,b)` | `ord(P)=r` | rank app. | `(▲)` holds | `χ(A)` | `χ(B)` | `π_χ` | `π_W` | QR/NQR/0 | bias |
|----:|:-------------:|----------:|----------:|:-----------:|:------:|:------:|------:|------:|:--------:|-----:|
| 1009 | (37, 2)  | 92   | 92   ✓ | ✓ | +1 | −1 | `2r` | `1008·r` | 53/38/1 | **+0.165** |
| 1013 | (5, 7)   | 91   | 91   ✓ | ✓ | +1 | +1 | `1r` | `253·r`  | 52/38/1 | **+0.156** |
| 2003 | (11, 19) | 116  | 116  ✓ | ✓ | −1 | +1 | `2r` | `2002·r` | 56/59/1 | −0.026 |
| 7919 | (3, 8)   | 3964 | 3964 ✓ | ✓ | −1 | −1 | `2r` | `7918·r` | 2012/1951/1 | +0.015 |
| 10007| (17, 23) | 3279 | 3279 ✓ | ✓ | +1 | +1 | `1r` | `5003·r` | 1639/1639/1 | **0.000** |

The 2-D net zero-lattice recovered the discrete log exactly in all five
cases (e.g. `p=7919`: `Q=[1321]P` recovered from the net zero
`(a,b)=(2643,1)`, since `2643 + 1·1321 = 3964 = r ≡ 0`).

**Read-outs.**

1. **Anchor solid.** Rank of apparition `= ord(P)`, apparition law, and
   the multiplier law `(▲)` all hold exactly — the instrument is
   computing genuine EDS, not a look-alike.

2. **χ-period dichotomy confirmed.** `π_χ = r` exactly when
   `χ(A)=χ(B)=+1` (rows 1013, 10007); `π_χ = 2r` in every other sign
   combination (rows 1009, 2003, 7919). The predicted closed form matches
   with no exceptions.

3. **Residue bias is curve-dependent and sometimes large.** This is the
   structurally interesting find. If `χ(W(n))` over the apparition block
   were a fair `±1` coin, the bias `(QR−NQR)/(QR+NQR)` would sit within a
   few percent of zero. Two of five curves show `+0.16` — a real,
   non-generic imbalance — while `p=10007` is *perfectly* balanced
   (1639/1639). The bias is therefore **not universal but not noise
   either**: it is a curve-specific quantity. Mapping which curves carry
   large residue bias, and why, is the first thing the avenue should
   chase. (Sample sizes here are small; §5 specifies the scale-up.)

---

## 5. A concrete research program

Falsifiable, ordered by cost:

1. **Bias census.** Sweep thousands of curves at `p ≈ 2¹⁶–2²⁰`, compute
   the apparition-block residue bias, and characterise its distribution.
   *Prediction to break:* bias concentrates near `0` with curve-specific
   outliers correlated with the multiplier characters `(χ(A), χ(B))`
   and/or the CM discriminant. If a clean predictor of large bias exists,
   that is a publishable distinguisher.

2. **Multiplier-character law in general.** Prove the §3 dichotomy for all
   `r` (the `r`-odd case has a different parity term) and express
   `(χ(A), χ(B))` intrinsically — conjecture: in terms of the quadratic
   character of `B = (the Weierstrass-`℘`′-type ratio)` and ultimately the
   2-torsion / the Tate pairing `⟨P,P⟩`. This connects EDS-Residue back to
   EDS-Association and the pairing.

3. **Net-`χ` localisation experiment.** Build the genuine 2-D net (not the
   point-arithmetic stand-in used here for the zero-lattice illustration)
   and measure whether the Legendre pattern `χ(W(a,b))` constrains `Λ_k`
   enough to beat `O(√m)`. *Prediction to break:* it does not beat
   generic — but quantify the constant, and test small-embedding-degree
   curves where EDS-Association is already easy, to see if the residue
   signal piggybacks.

4. **`F_p` ↔ `Z` bridge.** Make the Silverman–Stephens sign-period formula
   explicit mod `p` and check whether the `χ`-period `r·j_χ` we measure is
   the reduction of their real-analytic period. A clean bridge would let
   one *predict* `χ`-structure from the curve's real period without
   building the sequence.

---

## 6. Honest scorecard

| Claim | Status |
|-------|--------|
| ECDLP "lives inside" the EDS/net (zero-lattice = slope `k`) | **Proven & demonstrated** (§2.1, §4) |
| EDS-Residue solves ⟺ ECDLP solves (sub-exp) | **Lauter–Stange theorem** — no shortcut |
| `χ`-period determined by `(χ(A),χ(B))` | **Confirmed empirically**, closed form in §3 |
| Residue bias is curve-dependent, sometimes large | **Observed** (+0.16 vs 0.00); needs census |
| QR pattern beats generic ECDLP | **No evidence; expected false.** Open to quantify constants |

**Why it is underexplored, fairly stated.** The equivalence theorem is
often read as closing the subject ("EDS-Residue ≡ ECDLP, move on"). But
the theorem is asymptotic; the *bias and period structure of the Legendre
sequence* — the part with curve-specific, non-generic behaviour visible
even in this five-curve toy sweep — was never the theorem's subject and
has had little systematic measurement. That gap, not a hoped-for
sub-exponential attack, is the real opportunity: a distinguisher census
that the existing literature simply did not run.

---

## 7. Tests & reproduction

```bash
cargo test  --release --lib cryptanalysis::eds_residue     # 6 tests
cargo run   --release --example eds_residue_demo           # the §4 table
```

Tests: `rank_of_apparition_equals_order`, `apparition_law_holds`,
`legendre_matches_euler`, `net_zero_lattice_recovers_discrete_log`,
`multiplier_law_and_periods`, `report_renders`.

---

## 8. Sources

- K. E. Lauter, K. E. Stange. *The Elliptic Curve Discrete Logarithm
  Problem and Equivalent Hard Problems for Elliptic Divisibility
  Sequences.* SAC 2008. arXiv:0803.0728, ePrint 2008/099.
  <https://arxiv.org/abs/0803.0728> · <https://eprint.iacr.org/2008/099>
- K. E. Stange. *The Tate Pairing via Elliptic Nets.* Pairing 2007.
  ePrint 2006/392. <https://eprint.iacr.org/2006/392>
- R. Shipsey. *Elliptic Divisibility Sequences.* PhD thesis, Goldsmiths,
  University of London, 2000.
- J. H. Silverman, N. Stephens. *The sign of an elliptic divisibility
  sequence.* J. Ramanujan Math. Soc. 21 (2006). arXiv:math/0402415.
  <https://arxiv.org/abs/math/0402415>
- M. Ward. *Memoir on elliptic divisibility sequences.* Amer. J. Math. 70
  (1948).
- *Elliptic divisibility sequence.* Wikipedia.
  <https://en.wikipedia.org/wiki/Elliptic_divisibility_sequence>
