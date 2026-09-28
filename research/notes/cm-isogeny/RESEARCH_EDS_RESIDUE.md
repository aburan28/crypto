# Elliptic divisibility sequences, elliptic nets, and the EDS-Residue handle on the ECDLP

**Modules:** `src/cryptanalysis/eds_residue.rs`, `eds_tate.rs` (§5.6), `eds_net.rs` (§5.3b), `eds_mov.rs` (§5.10)
**Demos:** `examples/eds_residue_demo.rs`, `eds_census.rs`, `eds_localisation.rs`, `eds_bridge.rs`, `eds_tate_demo.rs`, `eds_net_demo.rs`, `eds_mov_demo.rs`
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

2. **Multiplier-character law. DONE.** Closed forms `(CF)` + structural
   identity `(BR)` (§5.5), and the sequence-free **Tate-pairing bridge**
   `χ(B) = χ(⟨P,P⟩_r)` confirmed in the nondegenerate regime (§5.6).

3. **χ-localisation: how many residue bits identify the discrete log?**
   **DONE: rank-1 (§5.3a) and rank-2 net (§5.3b, built without Stange's seeds).**

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

### 5.3b Genuine 2-D net — built without Stange's seeds, and the verdict

**Module:** `src/cryptanalysis/eds_net.rs` · **demo:** `eds_net_demo`.

Stange's mixed initial block (`W(2,1)`, `W(1,2)`, …; Props 6.3/6.4 of
arXiv:0710.1316) could not be fetched here, and the net recurrence alone
underdetermines them. So instead of guessing, the net is **derived** from a
rank-2 generalisation of the rank-1 coordinate relation
`ψ_{a+1}ψ_{a−1}/ψ_a² = x(P)−x(aP)`:

```
W(a+1,b)·W(a−1,b) / W(a,b)²  =  x(P) − x(aP+bQ)            (REL-P)
W(a,b+1)·W(a,b−1) / W(a,b)²  =  x(Q) − x(aP+bQ)            (REL-Q)
```

With gauge `W(1,0)=W(0,1)=W(1,1)=1`, axes `W(a,0)=ψ_a(P)`, `W(0,b)=ψ_b(Q)`,
and `x(aP+bQ)` from point arithmetic, these second-order recurrences fill the
grid — *no external seeds needed*. The result is then **validated**, not
assumed:

- the net recurrence `W(p+q)W(p−q)W(r)² = W(p+r)W(p−r)W(q)² − W(q+r)W(q−r)W(p)²`
  holds on every checked triple (`net_satisfies_recurrence`);
- the zero set is exactly `{(a,b) : aP+bQ = O}` (`net_zero_lattice_matches_point_arithmetic`);
- the axes reproduce the rank-1 EDS;
- **and it is genuinely rank-2**: on a full-7-torsion curve over `F_1009`
  with *independent* `P, Q` (`Q ∉ ⟨P⟩`), `(NET)` still holds and the
  zero-lattice is the 2-D sublattice `7ℤ×7ℤ`, not a line
  (`net_recurrence_holds_for_independent_p_q`). So the construction is *the*
  rank-2 net, not a degenerate artifact of `Q ∈ ⟨P⟩`.

So **§5.3b is unblocked**: a genuine canonical net is in hand, validated in
both the degenerate (ECDLP) and the truly 2-dimensional case.

**The verdict on χ-localisation.** For the ECDLP, `Q = [k]P`, so every point
`aP+bQ = [(a+bk) mod m]P` lives in `⟨P⟩`, and `(REL-P)` only ever consumes
`x(jP)` values — **the rank-2 net is a reparametrisation of the rank-1 EDS of
`P`**. Its `χ(W(a,b))` pattern therefore carries no localisation power beyond
§5.3a: information-tight (`~log₂ m` bits) but algorithmically inert
(`O(m)` extraction, no sub-`√m` advantage). The demo shows row `b=0` is
literally the rank-1 EDS χ-row, and confirms `aP+bQ=[(a+bk) mod m]P` for the
whole grid. A *non-degenerate* rank-2 net (independent `P,Q` in different
subgroups) does **not** arise in the ECDLP — so the 2-D net offers no opening
the 1-D handle didn't, consistent with Lauter–Stange. The EDS-Association ↔
Tate-pairing ↔ MOV/Frey–Rück route (§2.2, realised by
`cryptanalysis::mov_attack`) remains the only regime where this structure
collapses the DLP, exactly when the embedding degree is small.

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

### 5.6 The Tate-pairing bridge: `χ(B) = χ(⟨P,P⟩_r)` (confirmed)

**Module:** `src/cryptanalysis/eds_tate.rs` · **demo:** `eds_tate_demo`.

`(CF)`/`(BR)` still read terms *near index `r`*, so they need the `O(r)`
sequence. The sequence-free route goes through `B² = −W(r+1)/W(r−1)` and the
**self-Tate–Lichtenbaum pairing** `⟨P,P⟩_r`. For embedding degree 1
(`r ∣ p−1`) this pairing lives in `μ_r ⊂ F_p^*` and is computable by Miller's
algorithm in `O(log r)` — *without* the EDS. I implemented a full `F_p` Tate
pairing from scratch (Miller + final exponentiation, with the loop-ending
2-torsion / inverse-point cases that arise for even `r`) and **validated it
independently**: on every instance it lands in `μ_r` (`t^r=1`) and is
bilinear (`⟨P,[c]P⟩ = ⟨P,P⟩^c`). 84 instances over 6 primes, all valid.

The naive conjecture `χ(B) = χ(⟨P,P⟩_r)` holds only 62/84 — but the failures
are **structural, not random**. For `t ∈ μ_r`, `χ(t) = t^{(p−1)/2}` is forced
to `+1` whenever `v₂(r) < v₂(p−1)` (then `−1 ∉ μ_r`). Splitting by that
2-adic valuation gives a clean dichotomy:

| regime | meaning | result |
|--------|---------|--------|
| `v₂(r) = v₂(p−1)` (**nondegenerate**) | `−1 ∈ μ_r`, `χ(t)` can be `−1` | **`χ(B) = χ(⟨P,P⟩_r)`: 27/27** |
| `v₂(r) < v₂(p−1)` (forced) | `χ(t) = +1` carries no information | `χ(t)=+1`: 57/57 |

So the conjecture, *correctly stated*, is **confirmed**: whenever the
self-Tate pairing's quadratic character is non-degenerate
(`v₂(r) = v₂(p−1)`), it **equals** the EDS-Residue multiplier character
`χ(B)`. This is the bridge the program sought — it ties the EDS-Residue
χ-structure (§3, §5.5) directly to a **Tate pairing** (the object of
EDS-Association, §2.2), and in that regime `χ(B)` is computable in
`O(log r)` from the pairing without ever building the `O(r)` sequence.
(Tests: `tate_pairing_is_valid`,
`chi_b_equals_chi_self_tate_in_nondegenerate_regime`.)

### 5.7 The forced regime, resolved: the *unreduced* pairing character

The §5.6 bridge used the *reduced* pairing `t = f^{(p−1)/r} ∈ μ_r`. In the
forced regime `(p−1)/r` is even, so `χ(t) = χ(f)^{(p−1)/r} = +1` — the final
exponentiation destroys the quadratic character. The fix is to **not reduce**:
keep the unreduced Miller value `f = f_{r,P}(D_P)`. The Tate pairing fixes `f`
only up to `(F_p^*)^r`, but for **even `r`** every `r`-th power `z^r =
(z^{r/2})²` is a square, so `χ` is constant on those cosets — hence `χ(f)` is
a **well-defined, `S`-independent** bit. The result (test-verified, 84
even-`r` instances over 7 primes, all valid):

> **For every even `r`,  `χ(B) = χ(f_{r,P}(D_P))`** — the unreduced
> self-Tate-pairing character — *in both regimes*.

This subsumes §5.6: when `(p−1)/r` is odd (nondegenerate) `χ(f) = χ(t)`, so it
reduces to the earlier statement; when `(p−1)/r` is even (forced) `χ(t)` is
trivial but `χ(f)` still equals `χ(B)`. Two checks pin it down: `χ(f)` is
identical for two independent auxiliary points `S` (S-independence, as the
`r`-even coset argument predicts), and in the forced regime `χ(t) = +1`
always while `χ(f) = χ(B)` varies. So the EDS-Residue multiplier character is
the quadratic character of the (unreduced) self-Tate pairing for *all* even
`r` — and it stays `F_p`-computable; no exotic handle is needed.
(Test: `unreduced_self_tate_char_equals_chi_b_all_even_r`.)

**What is left.** Only odd `r` sits outside this pairing statement (there
`(F_p^*)^r ⊄` squares, so `χ(f)` is ill-defined) — but odd `r` already has
the clean closed form `χ(B)=χ(−1)χ(W(r+1))χ(W(r−1))` from §5.5. Between §5.5
(odd `r`) and §5.7 (even `r`), `χ(B)` is now pinned in every case.

### 5.8 Lift to embedding degree > 1 (the MOV regime), still in `F_p`

§5.6/§5.7 assumed embedding degree 1 (`r ∣ p−1`), so the pairing lived in
`F_p`. The genuine MOV/Frey–Rück regime has embedding degree `k > 1`: `r ∣
p^k−1` but `r ∤ p−1`, and the *reduced* Tate pairing lands in
`μ_r ⊂ F_{p^k}^*`. One might expect the bridge to need extension-field
arithmetic there. It does not. The **unreduced** self-Miller value
`f_{r,P}(D_P)` is computed from the `F_p`-rational Miller function evaluated
at `F_p` points, so **it stays in `F_p` no matter the embedding degree** —
only the final exponentiation would leave `F_p`. The result (test-verified on
≥20 embedding-degree-`>1` instances, even `r`, over 6 primes):

> **For every even `r`, any embedding degree, `χ(B) = χ(f_{r,P}(D_P))`**,
> with `χ(f)` `S`-independent — entirely within `F_p`.

(For even `r`, `gcd(r,p−1)` is even, so `(F_p^*)^r ⊆` squares and `χ(f)` is
well-defined on the divisor-ambiguity cosets regardless of `k`.) So the
EDS-Residue ↔ self-pairing identification is **not** a low-embedding-degree
accident: it is an `F_p` statement that holds across the whole MOV spectrum.
The only place the embedding degree matters is whether the *reduced* pairing
(hence an actual sub-exponential DLP transfer, MOV) is available — which is
the standard `k`-small criterion, orthogonal to the χ-bridge.
(Test: `unreduced_self_miller_char_equals_chi_b_any_embedding_degree`.)

The sharpest case is the **supersingular** family `y²=x³+x` over `p≡3
(mod 4)`: `#E=p+1`, embedding degree *exactly* 2 — the canonical MOV-weak /
pairing curves. The bridge `χ(B)=χ(f_{r,P})` holds there too, computed wholly
in `F_p`, S-independent, across 6 primes and several orders per curve. So the
identification survives precisely on the curves where the *reduced* pairing
would transfer the DLP into `F_{p²}` (Test:
`bridge_holds_on_supersingular_curves`).

### 5.9 Stange's Tate-pairing-via-net formula, cross-validated

The non-degenerate net (§5.3b) finally lets us close the loop with Stange's
original motivation — computing the Tate pairing as a *ratio of net values*.
For `P` of order `r` and an arbitrary `Q`, Stange's formula (embedding degree
1) is

```
τ_r(P,Q) = ( W(r+1,1)·W(1,0) / (W(r+1,0)·W(1,1)) )^{(p−1)/r}.
```

Crucially, when `P,Q` are **independent** (`Q ∉ ⟨P⟩`) the row `b=1` has *no*
zeros — `aP+Q = O` is impossible — so `W(r+1,1)` is reachable by the
`(REL-P)` fill even though the axis row `b=0` vanishes at `a=r`. On a full
7-torsion curve over `F_1009` (so `r=7 ∣ p−1`, independent order-7 `P,Q`):

> the net ratio `τ_net` **equals** the independent Miller-based Tate pairing
> `τ_miller` exactly, and both are a **nondegenerate** primitive `r`-th root
> of unity (`≠1`, in `μ_r`).

This is a **triple cross-validation** in one identity: it confirms (i) the
`(REL-P)/(REL-Q)` net is the canonical Stange net (a *wrong* gauge would not
satisfy the gauge-invariant Tate ratio), (ii) the from-scratch Miller pairing
of §5.6, and (iii) Stange's net-Tate formula itself. The EDS / elliptic-net
machinery built here is therefore mutually consistent end to end — the net
computes the pairing, the pairing's character is the EDS multiplier character
(§5.7–§5.8), and the multiplier governs the χ-structure (§3–§5.5).
(Test: `tate_via_net_matches_miller`.)

### 5.10 The payoff: a working MOV attack via the pairing

Everything above says the EDS/net handle is *inert* on generic curves — but
it **computes the pairing**, and on the right curves the pairing is not
inert. `eds_mov.rs` closes the arc with an end-to-end MOV/Frey–Rück break.

On the supersingular curve `y²=x³+x` over `p≡3 (mod 4)` (`#E=p+1`, embedding
degree 2), the distortion map `φ(x,y)=(−x, iy)` (`i²=−1 ∈ F_{p²}`) sends `P ∈
E(F_p)` to an independent point of `E(F_{p²})`, so the modified reduced Tate
pairing `t_r(P, φ(P)) = f_{r,P}(φ(P))^{(p²−1)/r} ∈ μ_r ⊂ F_{p²}^*` is
**non-degenerate**. Built from scratch here (minimal `F_{p²} = F_p[i]`
arithmetic + Miller evaluated at the `F_{p²}` point `φ(P)`) and validated:
the pairing lands in `μ_r` (`t^r=1`), is bilinear (`t([c]P,φ(P)) = t^c`), and
`≠1`. Bilinearity then transfers the ECDLP `Q=[k]P` into a DLP `α^k = β` in
the *small* cyclic group `μ_r`, solved generically:

```
E: y²=x³+x / F_1283  (supersingular, #E=1284), r=107 | p+1
P=(921,458) ord 107,  t_r(P,φ(P))=(171,824) ∈ μ_107  (nondegenerate)
Q=[7]P→k=7   Q=[42]P→k=42   Q=[99]P→k=99   Q=[106]P→k=106   ✓
```

This is the concrete edge of the whole investigation. The EDS-Residue /
elliptic-net handle gives **no** sub-`√m` attack on generic curves (§3–§5.9),
but the same machinery *is* the Tate pairing (§5.6, §5.9), and the pairing
collapses the ECDLP **exactly** on the small-embedding-degree curves — which
is precisely the Lauter–Stange "EDS Association ⇒ weak curve" boundary (§2.2)
made executable. (Tests: `distortion_point_is_on_curve_over_fp2`,
`pairing_is_nondegenerate_and_bilinear`, `mov_attack_recovers_discrete_log`.)

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
| Multiplier characters have closed forms `(CF)` + identity `Bʳ=−W(r+1)W(r−1)` | **Derived & test-verified (§5.5)** on 5/3 curves |
| Tate bridge `χ(B) = χ(⟨P,P⟩_r)` when `v₂(r)=v₂(p−1)` | **Confirmed (§5.6):** 27/27 nondeg; F_p Tate pairing built & validated (bilinear, μ_r), 84 instances |
| Forced regime `v₂(r)<v₂(p−1)`: `χ(B)` from the *unreduced* pairing | **Resolved (§5.7):** `χ(B)=χ(f_{r,P})` for all even `r`, both regimes, S-independent — test-verified |
| Bridge lifts to embedding degree > 1 (MOV regime) | **Yes (§5.8):** `χ(B)=χ(f_{r,P})` for `r∤p−1`, S-independent, *entirely in `F_p`* — incl. **supersingular** `y²=x³+x` (embedding degree 2), test-verified |
| Canonical 2-D net derivable without Stange's seeds | **Yes (§5.3b):** built via (REL-P)/(REL-Q), validated by (NET) + zero-lattice + axes |
| Stange's Tate-via-net formula reproduces the Miller pairing | **Yes (§5.9):** `τ_net = τ_miller` (nondegenerate, in `μ_r`) — triple cross-validation of net, pairing, formula |
| The pairing *does* break the ECDLP on weak curves | **Yes (§5.10):** end-to-end MOV attack on supersingular `y²=x³+x` via `F_{p²}` distortion pairing — recovers `k` |
| QR pattern (1-D or 2-D net) beats generic ECDLP | **No.** Info-tight but algorithmically inert (§5.3a); for `Q∈⟨P⟩` the 2-D net is a rank-1 reparametrisation (§5.3b) — no sub-`√m` advantage |

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
holds; the χ-structure is reduced to closed forms in the multiplier (§5.5);
and that multiplier character is **identified with a self-Tate pairing** —
`χ(B)=χ(⟨P,P⟩_r)` in the nondegenerate regime (§5.6), and more generally
`χ(B)=χ(f_{r,P})` (the unreduced pairing character) for every even `r`
including the forced regime (§5.7) and *every embedding degree* (§5.8, the
MOV regime, still entirely in `F_p`) — the concrete bridge from EDS-Residue to
EDS-Association the program set out to find; and the
canonical 2-D net, the last "blocked" item, was **derived from scratch**
without Stange's seeds (§5.3b) and shown to be a rank-1 reparametrisation in
the ECDLP case — no new opening. Every item of the program is now settled:
the EDS/elliptic-net handle on the ECDLP is real, rigidly structured, and
**information-tight but algorithmically inert** — no sub-`√m` attack, exactly
as Lauter–Stange's equivalence predicts, now demonstrated end to end.

---

## 7. Tests & reproduction

```bash
cargo test  --release --lib cryptanalysis::eds_residue     # 20 tests
cargo test  --release --lib cryptanalysis::eds_tate        #  5 tests (§5.6–§5.8)
cargo test  --release --lib cryptanalysis::eds_net         #  4 tests (§5.3b, §5.9)
cargo test  --release --lib cryptanalysis::eds_mov         #  3 tests (§5.10)
cargo run   --release --example eds_residue_demo           # the §4 table
cargo run   --release --example eds_census                 # the §4.5 census
cargo run   --release --example eds_localisation           # the §5.3a sweep
cargo run   --release --example eds_bridge                 # the §5.4 bridge
cargo run   --release --example eds_tate_demo              # the §5.6 bridge
cargo run   --release --example eds_net_demo               # the §5.3b net
cargo run   --release --example eds_mov_demo               # the §5.10 MOV attack
```

`eds_mov` tests: `distortion_point_is_on_curve_over_fp2`,
`pairing_is_nondegenerate_and_bilinear`, `mov_attack_recovers_discrete_log`.

`eds_tate` tests: `tate_pairing_is_valid`,
`chi_b_equals_chi_self_tate_in_nondegenerate_regime`,
`unreduced_self_tate_char_equals_chi_b_all_even_r`,
`unreduced_self_miller_char_equals_chi_b_any_embedding_degree`,
`bridge_holds_on_supersingular_curves`.

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
