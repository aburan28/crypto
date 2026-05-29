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

3. **Residue bias looks curve-dependent at small `m`** — two curves show
   `+0.16`, while `p=10007` is *perfectly* balanced (1639/1639). At this
   sample size one cannot tell a structural distinguisher from finite-size
   noise. **The census in §4.5 resolves this**, and the answer is sharper
   than "curve-specific quantity": the balanced curve is exactly the
   `(χA,χB)=(+,+)` case, and there is an *exact symmetry law* behind it.

---

## 4.5 Bias census + the reflection-symmetry law (step 1, completed)

`cargo run --release --example eds_census` sweeps thousands of curves over
five primes — three with `p ≡ 3 (mod 4)` and two with `p ≡ 1 (mod 4)` — and
aggregates the apparition-block residue bias by multiplier class. The result
is a **clean, decisive answer** with a sharp falsifiable prediction that
checks out.

| `p` | class | curves | mean bias | std(bias) | mean \|bias\| | max \|bias\| | heavy tail |
|----:|:-----:|-------:|----------:|----------:|--------------:|-------------:|-----------:|
| 4 099 | ≡3 | 1 568 | +0.0013 | 0.0361 | 0.0204 | 0.316 | 7.7 % |
| 10 007 | ≡3 | 1 564 | +0.0012 | 0.0215 | 0.0130 | 0.150 | 8.0 % |
| 100 003 | ≡3 | 760 | +0.0001 | 0.0066 | 0.0039 | 0.045 | 8.3 % |
| 4 093 | ≡1 | 1 557 | +0.0011 | 0.0320 | 0.0189 | 0.302 | 2.7 % |
| 10 009 | ≡1 | 1 553 | +0.0005 | 0.0203 | 0.0124 | 0.174 | 6.1 % |

Mean \|bias\| per multiplier class — note the **balanced class flips with
`p mod 4`**:

| class `(χA,χB)` | 4099 (≡3) | 10007 (≡3) | 100003 (≡3) | 4093 (≡1) | 10009 (≡1) |
|:---------------:|----------:|-----------:|------------:|----------:|-----------:|
| `(+,+)` | **0.0000** | **0.0000** | **0.0000** | 0.0300 | 0.0191 |
| `(+,−)` | 0.0195 | 0.0123 | 0.0040 | 0.0188 | 0.0123 |
| `(−,+)` | 0.0335 | 0.0214 | 0.0066 | **0.0000** | **0.0000** |
| `(−,−)` | 0.0210 | 0.0128 | 0.0034 | 0.0212 | 0.0138 |

**Two facts jump out.**

(a) **The bias decays like `1/√m`.** As the order grows (`m ≈ 150 → 600 →
3 800`), `max|bias|` falls `0.316 → 0.150 → 0.045` and `std(bias)` falls
`0.036 → 0.021 → 0.0066` — the exact scaling of the sample standard
deviation of a fair `±1` coin over `m` draws. The heavy-tail fraction sits
at a constant `≈ 8 %`, i.e. ordinary Gaussian fluctuation (no structured
excess of outliers). **So there is no constant residue-bias distinguisher**
— the `+0.16` values at §4 were small-`m` noise, and they vanish at
cryptographic scale. This is the negative result the avenue most needed,
and it is now nailed down rather than assumed.

(b) **But the `(+,+)` class has bias *exactly* 0 at every scale, and every
single heaviest-bias curve is `(−,+)`.** That is not noise — it is an exact
symmetry. The EDS extends to negative indices by `W(−n) = −W(n)`, so the
multiplier law `(▲)` gives the **reflection identity**

```
W(r − n) = −A · B^{−n} · W(n)      ⇒
χ(W(n))·χ(W(r−n)) = χ(−1) · χ(A) · χ(B)ⁿ.                          (◆)
```

For the block to be *forced* balanced, the pair product `(◆)` must equal
`−1` for every `n`, which (since `χ(B)ⁿ` must be constant in `n`) happens
**iff `χ(B) = +1` and `χ(A) = −χ(−1)`**. Concretely:

- **`p ≡ 3 (mod 4)`** (`χ(−1) = −1`): balanced class is `(χA,χB) = (+,+)`;
  the *reinforcing* class (pair product `+1`, largest fluctuations) is
  `(−,+)`. Every census outlier at `p≡3` is `(−,+)`. ✓
- **`p ≡ 1 (mod 4)`** (`χ(−1) = +1`): the balanced class **flips to
  `(−,+)`**, and `(+,+)` becomes the reinforcing class. The census confirms
  this exactly — at `p = 4093, 10009` the `(−,+)` column is `0.0000` and
  every top-bias curve is `(+,+)`. ✓
- **`χ(B) = −1`:** the sign alternates with `n`; no global cancellation,
  intermediate behaviour (`(+,−)` / `(−,−)` rows).

(The `0` is exact up to the lone `n = r/2` fixed point when `r` is even.)
This `χ(−1)`-dependent flip is the kind of *predicted-then-confirmed* result
that distinguishes a real law from a fitted curiosity. Identity `(◆)` is
verified term-by-term in the test suite for **both** regimes
(`reflection_symmetry_law` at `p≡3`, `reflection_law_holds_for_p_eq_1_mod_4`
at `p≡1`), and the balanced-class consequence in `plus_plus_class_is_balanced`
and `balanced_class_flips_for_p_eq_1_mod_4`.

**Upshot.** The residue "bias" is *not* an exploitable distinguisher: it is
finite-size sampling noise (`∝ 1/√m`), modulated by an exact reflection
symmetry whose balanced class is `{χ(B)=+1, χ(A)=−χ(−1)}`. The Legendre
sequence is, at cryptographic scale, statistically generic — but its
*structure* (period dichotomy §3, reflection law `(◆)`) is rigid and fully
predictable from `(χ(A), χ(B))` and `p mod 4`. The handle is real; the
*bias-based* attack on it is closed.

---

## 5. A concrete research program

Falsifiable, ordered by cost:

1. ~~**Bias census.**~~ **DONE (§4.5), both residue classes.** Outcome: the
   bias is finite-size noise (`∝ 1/√m`, ~3–8 % Gaussian tail), with **no**
   constant distinguisher, but governed by the exact reflection law `(◆)`
   whose balanced class `{χ(B)=+1, χ(A)=−χ(−1)}` flips with `p mod 4` —
   predicted and confirmed (`(+,+)` balanced at `p≡3`, `(−,+)` at `p≡1`).
   The hoped-for "clean predictor of large bias" does not exist as an
   *advantage*; it exists as a *symmetry*. Remaining only: push to
   `p ≈ 2³²` to confirm the `1/√m` extrapolation (the u64 path supports it;
   cost-bounded by `min_order`/`cap`).

2. **Multiplier-character law. DONE for the closed form (§5.5); the
   sequence-free / Tate-pairing step remains a stated lead.**

3. **χ-localisation: how many residue bits identify the discrete log?**
   **DONE in rank-1 (§5.3a); rank-2 net deferred (§5.3b).**

### 5.3a Rank-1 χ-localisation via the decimation identity

The genuine 2-D net is not required to ask the sharp EDS-Residue question.
The **decimation identity** (Ward / Shipsey)

```
ψ_{κn}(P) = ψ_n([κ]P) · ψ_κ(P)^{n²}                                  (DEC)
```

gives, on Legendre symbols (using `χ(·)^{n²} = χ(·)ⁿ`),

```
χ(ψ_n(Q)) = χ(ψ_{kn}(P)) · χ(ψ_k(P))ⁿ      for Q = [k]P.            (DEC-χ)
```

The left side is computable from `Q`'s coordinates **alone** (no knowledge
of `k`); every `n` is therefore one EDS-Residue bit constraining `k`
against the public table `T[j] = χ(ψ_j(P))`. `localisation_sweep` measures
the minimal window of `n`'s that pins `k`. `cargo run --release --example
eds_localisation`:

| `p` | `p mod 4` | `ord(P)` | tested `k` | pinned to `±k` | median window | max window | sign resolved |
|----:|:---------:|---------:|-----------:|:--------------:|:-------------:|:----------:|:-------------:|
| 4099 | 3 | 64 | 61 | 61/61 | 6 | 10 | **61/61** |
| 2003 | 3 | 116 | 113 | 113/113 | 8 | 11 | **113/113** |
| 1019 | 3 | 524 | 200 | 200/200 | 10 | 23 | **200/200** |
| 1009 | 1 | 92 | 89 | 89/89 | 7 | 11 | **0/89** |
| 2017 | 1 | 694 | 200 | 200/200 | 10 | 17 | **0/200** |
| 4093 | 1 | 2029 | 200 | 200/200 | 11 | 23 | **0/200** |

Three clean read-outs:

- **The residues always pin `k`** (up to the unavoidable `±k`): `100 %`
  on every curve. `DEC-χ` verified for every `k` along the way.
- **`log₂ m` scaling.** Median window `6 → 7 → 8 → 10 → 10 → 11` as
  `m: 64 → 92 → 116 → 524 → 694 → 2029` tracks `log₂ m` (≈ 6 … 11). Each
  residue contributes ≈ one bit, and `k` carries `log₂ m` bits — the signal
  is **information-theoretically tight**, no redundancy, no obstruction.
- **A third `p mod 4` dichotomy.** The sign of `k` is resolved iff
  `χ(−1) = −1`, i.e. `p ≡ 3 (mod 4)` (`100 %` resolved), and is *never*
  resolved at `p ≡ 1` (`0 %`). Reason: `ψ_n(−Q) = ±ψ_n(Q)` with the sign
  `= χ(−1)` on even `n`, so the `k ↔ m−k` pair is distinguishable exactly
  when `χ(−1) = −1`.

**But this is not a computational shortcut.** Using the bits still requires
the candidate scan — `O(m)` survivors per window, `O(m·log m)` total —
which is *worse* than the generic `O(√m)`. So the EDS-Residue signal is
**information-rich yet algorithmically inert**: a textbook illustration of
*why* Lauter–Stange's equivalence holds — the bits are there and tight, but
extracting `k` from them costs a full search. (Tests:
`decimation_identity_holds`, `localisation_pins_the_true_k`,
`localisation_sweep_sign_dichotomy_and_log_window`.)

### 5.3b Genuine 2-D net (deferred)

The remaining piece — does the Legendre pattern of the *canonical 2-D net*
`χ(W(a,b))` localise `Λ_k` cheaper than `O(√m)`? — needs the net's mixed
initial block (`W(2,1)`, `W(1,2)`, `W(2,−1)`; Props 6.3/6.4 of Stange,
arXiv:0710.1316). The fundamental recurrence
`W(p+q)W(p−q)W(r)² = W(p+r)W(p−r)W(q)² − W(q+r)W(q−r)W(p)²` is in hand and
the axes are the 1-D EDS this module validates, but the recurrence alone
underdetermines the mixed seeds, and they could not be retrieved here
(arXiv/ePrint returned HTTP 403). It is left unimplemented rather than
guessed — a wrong-but-zero-preserving gauge would pass a zero check yet
report a meaningless `χ`-pattern. **Certification plan once the seeds are
in hand:** (i) match the recurrence on both axes against the 1-D EDS (pins
the gauge), (ii) reproduce every net zero against `[a]P+[b]Q=O`, (iii)
match the mixed seeds to their closed forms. Note the EDS-Association ↔
Tate-pairing ↔ MOV/Frey–Rück link (§2.2) is *already* realised in this
crate by `cryptanalysis::mov_attack` — the regime where the net signal
provably collapses to an easy DLP.

4. **`F_p` ↔ `Z` bridge. DONE (§5.4).** Answer: the two periods are
   *orthogonal*. The archimedean sign is one fixed aperiodic object; the
   `F_p` χ-period hops between `r` and `2r` with `p` per the §3 law. No
   reduction relates them.

### 5.4 The archimedean sign (Silverman–Stephens) vs the F_p χ-period

A natural hope is that the `F_p` χ-period `r·j_χ` (§3) is the mod-`p`
shadow of the *archimedean* sign-period that Silverman–Stephens (2006)
attach to an integer EDS. It is not — and the bridge experiment shows why,
fully validated.

`eds_integer` builds the genuine **integer** EDS of curve `37a`
(`y²+y = x³−x`), point `(0,0)`, from the seeds `W(2),W(3),W(4) = 1,−1,1`
via the duplication formulas over `Z`. It reproduces **OEIS A006769**
exactly (`0,1,1,−1,1,2,−1,−3,−5,7,−4,−23,29,…`), test-checked to `n=25`.

- **Archimedean side.** `37a` has discriminant `> 0` (two real components)
  and an irrational rotation number, so by Silverman–Stephens its sign
  sequence is **aperiodic**. Measured: `++-++---+--+++--+---++-+++--+--+++-++---…`,
  no period `≤ 80` over 220 terms (`integer_eds_signs_aperiodic_for_37a`).
  This is a single fixed real-analytic object — the placement of `nP` on
  the real components, governed by the elliptic logarithm of `P`.

- **Arithmetic side.** Reducing the *same* integer sequence mod `p` and
  reading its F_p structure (`reduce_and_analyze`): the rank of apparition
  equals `ord(P mod p)` — verified against independent point arithmetic on
  the short form `Y²=x³−x+1/4`, `P=(0,1/2)`, for `p ∈ {7,11,13,23,29}` — and
  the χ-period is `r·j_χ` with `j_χ ∈ {1,2}` set by the multiplier
  characters (§3), **jumping with `p`**:

  | `p` | `ord(P mod p)` | `(χA,χB)` | χ-period |
  |----:|---------------:|:---------:|:--------:|
  | 7 | 9 | (−,+) | `2r=18` |
  | 13 | 16 | (+,−) | `2r=32` |
  | 23 | 11 | (+,+) | `1r=11` |
  | 29 | 12 | (+,+) | `1r=12` |
  | 41 | 51 | (+,+) | `1r=51` |
  | 43 | 14 | (+,−) | `2r=28` |

**Conclusion.** The archimedean sign-period (real, fixed, aperiodic) and the
`F_p` χ-period (arithmetic, `p`-dependent, `r` or `2r`) are different
invariants of the same point — one set by the *real* elliptic logarithm,
the other by the *quadratic character mod `p`* of the multiplier. There is
no naive reduction from one to the other, so the χ-structure cannot be read
off the curve's real period; it must be computed mod `p`. (Tests:
`integer_eds_matches_oeis_a006769`, `integer_eds_signs_aperiodic_for_37a`,
`bridge_reduction_matches_group_order_and_chi_law`.)

### 5.5 Intrinsic multiplier characters (program item 2)

The §3 χ-period and the §4.5 reflection law are both driven by the two sign
bits `(χ(A), χ(B))` of the shift multiplier `W(n+r) = A·Bⁿ·W(n)`. These bits
are not opaque. From the extraction `B = W(r+2)/(W(2)·W(r+1))` and
`A = W(r+1)/B = W(2)·W(r+1)²/W(r+2)`, and because `W(r+1)²` is a square:

```
χ(A) = χ( W(2)·W(r+2) ),     χ(B) = χ( W(2)·W(r+1)·W(r+2) ).          (CF)
```

`(CF)` is verified against the inversion-extracted `A, B` on 5 curves
(`multiplier_chars_closed_form_matches_extraction`). There is also a
**structural identity**, derived from Ward's relation at the rank of
apparition (`ψ_r(P)=0`) plus the periodicity `(▲)`:

```
Bʳ = − W(r+1)·W(r−1),     together with   A·B = W(r+1).               (BR)
```

`(BR)` is test-verified exactly mod `p` on 3 curves
(`structural_identity_b_pow_r`), sign and all. For **odd `r`**, `Bʳ` carries
the same character as `B` (since `B²` is a square), so `(BR)` collapses to
the near-intrinsic form

```
χ(B) = χ(−1)·χ(W(r+1))·χ(W(r−1)),     χ(A) = χ(−1)·χ(W(r−1))    (r odd),
```

also test-verified. These turn the entire χ-structure (period §3, balanced
class §4.5) into closed forms in three sequence values.

**What remains (the genuine lead).** `(CF)`/`(BR)` still read three terms
*near index `r`*, so they need the `O(r)` sequence. The sequence-free goal
runs through `B² = −W(r+1)/W(r−1)`: this ratio is conjecturally a
**self-pairing** value (Tate / Frey–Rück `⟨P,P⟩_r`), which for embedding
degree 1 (`r ∣ p−1`) lives in `μ_r ⊂ F_p^*` and is computable by Miller in
`O(log r)` — *without* the EDS. If `χ(B) = χ(⟨P,P⟩_r)` holds, the χ-period
becomes predictable from the pairing alone, **tying EDS-Residue directly to
EDS-Association** (the Lauter–Stange problem that *is* the Tate pairing,
§2.2). Testing this needs a correct `F_p` self-Tate-pairing; the crate's
only pairing is a k=2 pseudo-Weil routine (`mov_attack`), so it is left as
the next implementation step, with `(BR)` as the bridge to it.

---

## 6. Honest scorecard

| Claim | Status |
|-------|--------|
| ECDLP "lives inside" the EDS/net (zero-lattice = slope `k`) | **Proven & demonstrated** (§2.1, §4) |
| EDS-Residue solves ⟺ ECDLP solves (sub-exp) | **Lauter–Stange theorem** — no shortcut |
| `χ`-period determined by `(χ(A),χ(B))` | **Confirmed empirically**, closed form in §3 |
| Residue bias is an exploitable distinguisher | **Refuted (§4.5):** bias `∝ 1/√m`, ~8 % Gaussian tail, no constant advantage |
| Reflection law `(◆)`; balanced class `{χB=+1, χA=−χ(−1)}` flips with `p mod 4` | **Predicted & confirmed**, test-verified both regimes (§4.5), 5 primes |
| EDS residues pin `k` (up to `±`) in `~log₂ m` bits | **Measured (§5.3a):** 100 % of `k`, 6 primes; info-theoretically tight |
| Sign of `k` resolved iff `p ≡ 3 (mod 4)` | **Predicted & confirmed** (§5.3a): 100 % at `p≡3`, 0 % at `p≡1` |
| `F_p` χ-period = reduction of the archimedean sign-period | **Refuted (§5.4):** orthogonal invariants; integer EDS = A006769 (validated), signs aperiodic, χ-period hops `r`/`2r` with `p` |
| Multiplier characters have closed forms `(CF)` + identity `Bʳ=−W(r+1)W(r−1)` | **Derived & test-verified (§5.5)** on 5/3 curves; sequence-free Tate-pairing step open |
| QR pattern beats generic ECDLP | **No.** Residues are info-tight but algorithmically inert (`O(m log m)` scan ≫ `√m`, §5.3a). Canonical 2-D net (§5.3b) blocked on Stange's mixed seeds (arXiv 403 here) |

**Why it is underexplored, fairly stated.** The equivalence theorem is
often read as closing the subject ("EDS-Residue ≡ ECDLP, move on"). The
theorem is asymptotic, so it never spoke to the *finite statistics* of the
Legendre sequence — and that is the gap this note actually fills. The
honest result is three-sided: (i) the bias-distinguisher hope is
**refuted** (§4.5) — at cryptographic `m` the Legendre sequence is
statistically generic; (ii) the *structure* is rigid and was un-catalogued
— the period dichotomy (§3), the reflection law `(◆)` (§4.5), and the
`log₂ m`-tight localisation with its sign dichotomy (§5.3a) are exact,
test-verified, and predicted from `(χ(A),χ(B))` and `p mod 4`; and (iii)
the residues are **information-tight but algorithmically inert** — they
determine `k` in `~log₂ m` bits yet give no sub-`√m` algorithm, which is
the cleanest concrete illustration of *why* the Lauter–Stange equivalence
holds; and the χ-structure is now reduced to closed forms in the multiplier
(§5.5). Two related leads remain, both routing through a pairing this
sandbox can't yet supply: §5.3b (does the *canonical 2-D net's* `χ`-pattern
localise `Λ_k` cheaper than generic? — needs Stange's mixed seeds) and §5.5
(is `χ(B)` the character of the self-Tate-pairing `⟨P,P⟩_r`? — needs an
`F_p` Tate pairing). Both would tie EDS-Residue to EDS-Association/the Tate
pairing, which is the natural next thread.

---

## 7. Tests & reproduction

```bash
cargo test  --release --lib cryptanalysis::eds_residue     # 20 tests
cargo run   --release --example eds_residue_demo           # the §4 table
cargo run   --release --example eds_census                 # the §4.5 census
cargo run   --release --example eds_localisation           # the §5.3a sweep
cargo run   --release --example eds_bridge                 # the §5.4 bridge
```

Tests: `rank_of_apparition_equals_order`, `apparition_law_holds`,
`legendre_matches_euler`, `net_zero_lattice_recovers_discrete_log`,
`multiplier_law_and_periods`, `reflection_symmetry_law`,
`reflection_law_holds_for_p_eq_1_mod_4`, `plus_plus_class_is_balanced`,
`balanced_class_flips_for_p_eq_1_mod_4`, `sqrt_u64_roundtrips_both_residues`,
`decimation_identity_holds`, `localisation_pins_the_true_k`,
`localisation_sweep_sign_dichotomy_and_log_window`,
`integer_eds_matches_oeis_a006769`, `integer_eds_signs_aperiodic_for_37a`,
`bridge_reduction_matches_group_order_and_chi_law`,
`multiplier_chars_closed_form_matches_extraction`, `structural_identity_b_pow_r`,
`census_runs_and_is_sane`, `report_renders`.

---

## 8. Sources

- K. E. Lauter, K. E. Stange. *The Elliptic Curve Discrete Logarithm
  Problem and Equivalent Hard Problems for Elliptic Divisibility
  Sequences.* SAC 2008. arXiv:0803.0728, ePrint 2008/099.
  <https://arxiv.org/abs/0803.0728> · <https://eprint.iacr.org/2008/099>
- K. E. Stange. *The Tate Pairing via Elliptic Nets.* Pairing 2007.
  ePrint 2006/392. <https://eprint.iacr.org/2006/392>
- K. E. Stange. *Elliptic Nets and Elliptic Curves.* Algebra & Number
  Theory 5 (2011). arXiv:0710.1316. <https://arxiv.org/abs/0710.1316>
  (rank-2 net polynomials and initial values: Props 3.8, 6.3, 6.4 — needed
  for §5.3.)
- R. Shipsey. *Elliptic Divisibility Sequences.* PhD thesis, Goldsmiths,
  University of London, 2000.
- J. H. Silverman, N. Stephens. *The sign of an elliptic divisibility
  sequence.* J. Ramanujan Math. Soc. 21 (2006). arXiv:math/0402415.
  <https://arxiv.org/abs/math/0402415>
- M. Ward. *Memoir on elliptic divisibility sequences.* Amer. J. Math. 70
  (1948).
- *Elliptic divisibility sequence.* Wikipedia.
  <https://en.wikipedia.org/wiki/Elliptic_divisibility_sequence>
