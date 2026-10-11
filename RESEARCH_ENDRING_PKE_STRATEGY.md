# Research Note: Public-Key Encryption from the Endomorphism Ring Problem

**Date.** 2026-10-10. **Status.** Literature survey plus a design proposal. No code,
no measurements; every parameter figure below is a derived estimate, not a
measurement. **Goal addressed.** Devise a strategy for a *novel* public-key
encryption (PKE/KEM) scheme whose trapdoor is the endomorphism ring End(E) of
the public curve, building on primitives that exist today (Deuring
correspondence, Kani/higher-dimensional isogenies, lollipop attacks).
**Conductor.** `conductor check` could not reach a control plane in this
session (connection refused on localhost:8080); the note is filed under the
hook's warn-through path.

Sections 1-3 are established results with citations. Section 4 is the
proposal; the algebra there is ours and is marked as derived. Section 5 is
the cryptanalysis checklist and Section 6 the experiment protocol. Section 7
lists requirement status.

---

## 1. What the problem actually is

### 1.1 Hardness of EndRing (what the trapdoor is worth)

| Statement | Source | Status |
|---|---|---|
| EndRing ⇔ ℓ-IsogenyPath ⇔ MaxOrder (poly-time, GRH) | Wesolowski, FOCS 2021, ePrint 2021/919 | proved |
| EndRing ⇔ OneEnd ⇔ Isogeny (arbitrary degree), unconditional; worst-case ⇒ average-case | Page–Wesolowski EC 2024 (2023/1399); Herlédan Le Merdy–Wesolowski TCC 2025 (2025/271) | proved |
| Best classical: Õ(p^{1/2}) low memory; heuristic p^{1/3+o(1)} time *and* memory | EHLMP 2018; Wesolowski 2026/1486 (July 2026); Mamah 2026/1821 (concrete: p^{1/3} does not beat p^{1/2} at practical memory except maybe level I) | heuristic |
| Best quantum: Õ(p^{1/4}) (Grover / Biasse–Jao–Sankar); nothing subexponential for unstructured curves | SQIsign v2 spec §10.2.1; Mamah 2026/1821 | state of the art |
| One known non-scalar endomorphism α ⇒ classical |disc Z[α]|^{1/4}, quantum subexponential | Herlédan Le Merdy–Wesolowski CiC 2025 (2023/1448) | proved (GRH) |
| Any full-rank suborder, or an isogeny representation from a known-End curve ⇒ End in quantum poly time (IsERP) | Chen–Imran–Ivanyos–Kutas–Leroux–Petit AC 2023 (2023/779); Chen–Petit EC 2025 | proved |

Consequences for a design: the public key may be a curve and nothing else that
is "End-shaped". Publishing an endomorphism, an orientation, a suborder, or an
evaluable isogeny from a known-End curve collapses the hardness to
|disc|^{1/4} / subexponential / quantum-polynomial respectively. Also budget
p ≈ 2^{3λ} rather than 2^{2λ} if the 2026 p^{1/3} line is taken at face value;
the SQIsign round-3 primes were reportedly changed for this reason (2026/2221).

### 1.2 The structural obstruction

Given End(E_A) only, computing an isogeny E_A → E_C for a fresh E_C is
poly-time *equivalent* to computing End(E_C) (Page–Wesolowski Thm 8.5 and
Prop. 8.4; unconditional in 2025/271). So the owner of End(E_A) cannot invert
"send me a random curve isogenous to E_A". A PKE must put a **hint** in the
ciphertext, and the whole design question is: which hint is (a) useless to
everyone, (b) decodable with End(E_A), (c) cheap.

### 1.3 The only two known mechanisms for isogeny PKE

1. **Commutative group actions** (CSIDH, SCALLOP, PEGASIS/qt-Pegasis, CORAL).
   The orientation must be public for the sender to act, so vectorization is
   quantum-subexponential (Kuperberg); 2026/2215 puts NIST-1 at a 4096-bit
   base field. OSIDH tried a partially hidden orientation and was broken
   (Dartois–De Feo PKC 2022). Not EndRing-hard in the sense the user asked for.
2. **Torsion-point hints** (SIDH family: FESTA, QFESTA, POKÉ, PIKE, INKE,
   LIT-SiGamal, SILBE, Séta). Hints are torsion images hidden by a subgroup
   Γ ≤ GL₂(Z/N) (De Feo–Fouotsa–Panny EC 2024, 2024/459): Γ trivial or
   unipotent ⇒ polynomial for anyone (Robert 2022: N² > d suffices); scalar Γ
   (M-SIDH) and diagonal Γ (FESTA's CIST) ⇒ only exponential attacks known on
   generic curves.

Open-problem status: Castryck–De Feo–Galbraith–Kutas–Reijnders–Wesolowski,
"The isogeny problems" (ePrint 2026/1431), Problem 7 (De Feo): *build
encryption/KEM/IBE whose security reduces to a small variant of EndRing.*
Open as of Oct 2026. Problem 6 (Castryck): classical End transfer along an
arbitrary isogeny representation, also open.

## 2. Prior art that is "End(E) as trapdoor"

| Scheme | Trapdoor | Ciphertext hint | Why it is where it is |
|---|---|---|---|
| **Séta** (AC 2021, 2019/1291) | a curve built around a secret endomorphism θ (Petit's torsion attack as decryptor); needed N > D² | *exact* torsion images of a known-degree isogeny | broken 2022: Robert needs no θ once N² > d. Keygen 10.4 h, decryption 10.7 min in 2021. |
| **M-SIDH / MD-SIDH** (EC 2023, 2023/013) | key exchange, no trapdoor | scalar-masked images, N with ≥ 2λ prime factors | Weil pairing leaks λ² so only 2^{ω(N)} masks; needs E_0 of *unknown* End because any known endomorphism gives the lollipop ψ = φθφ̂ with the mask cancelled; 5911-bit p at level 1 |
| **Castryck–Vercauteren** (AC 2023, 2023/1433) | attack, not a scheme | generalised lollipop ψ = φ'∘ω∘φ̂, Lemma 3: the mask cancels iff the matrix of ω on E_0[N] lies in the centraliser of Γ; needs N > d√w | scalar masks: any ω works; F_p-rational E_0 is broken with unknown End (Frobenius). Diagonal masks (FESTA): only if P or Q is an *eigenvector* of ω. §5.3: a rigged eigenbasis is a *backdoor*, "near impossible to verify" when End(E_0) is unknown; they recommend hashing the basis. |
| **SILBE** (SAC 2024, 2024/400) | pk = any curve E_A; sk = φ_A: E_0 → E_A from F_p-rational E_0 (⇔ End(E_A)); decrypt with the Frobenius lollipop in dimension 4 | scalar-masked images (the mask *is* the message) | OW-PCA ⇔ "SIP + masked torsion"; first UPKE not from group actions; **p = 13 013 bits** at λ=128 because the scalar mask forces ω(N) ≥ 881 primes. Authors: "not yet practical". |
| **UPKE from FESTA** (AC 2026, 2026/1014) | dim-4 FESTA trapdoor; End of two curves used only for key update ("ConstrainedIsogeny") | diagonal-masked images | pk 0.4 KiB, ct 2.13 KiB; orders of magnitude faster than SILBE; assumption CIST² + IND-KLPT |
| FESTA / QFESTA / POKÉ / PIKE | secret isogeny from a *known-End* E_0 plus masks | diagonal (FESTA) or scalar (POKÉ) masks | key recovery is a torsion-hint problem, not EndRing. POKÉ: 431-bit p, pk 324 B, ct 648 B, ≈100 ms Sage. FESTA: 1292-bit p, pk 561 B, ct 1122 B. |
| pSIDH (Leroux 2022) | End(E_A) via suborder representation | the suborder itself | quantum-polynomial break (2023/779, Chen–Petit 2025) |

Reading of the table: the lollipop is the *only* published mechanism by which
End(E_A) decodes a hint that nobody else can decode, and its only published
instantiation (SILBE) is parametrically dead because it uses **scalar** masks,
whose security is 2^{ω(N)}. Diagonal masks have 2^{Ω(log N)} security but, by
CV Lemma 3, are lollipop-decodable only from an eigenbasis. CV treated the
eigenbasis as a backdoor to be eliminated. Nobody has used it as the
legitimate receiver's trapdoor. That is the gap this note proposes to fill.

## 3. The decoding lemma (why the design space is what it is)

Let φ: E_A → E_C have degree d, let (P,Q) be a basis of E_A[N], and let the
ciphertext reveal Λ = φ∘D on E_A[N], where D ∈ Γ ≤ GL₂(Z/N) is the sender's
mask written in the basis (P,Q). For θ ∈ End(E_A) with matrix M on E_A[N],
the receiver can compute

  Λ ∘ M ∘ (d·Λ⁻¹) = φ D M D⁻¹ φ̂   on E_C[N].

This equals the lollipop ψ_θ := φθφ̂ exactly iff D M = M D. (This is CV
Lemma 3 with σ₀ = id; the derivation above is a one-line restatement.) Hence:

* Γ = scalars: every θ works (M-SIDH, SILBE), but security is 2^{ω(N)}.
* Γ = diagonal in basis (P,Q): works iff M is diagonal, i.e. (P,Q) is an
  eigenbasis of θ. The centraliser of a non-scalar θ mod N is exactly the
  diagonal torus, so **diagonal masks are the largest mask set compatible
  with a single secret endomorphism**; two non-commuting endomorphisms would
  force scalars again.
* Once ψ_θ is known on E_C[N] with N² > deg ψ_θ = d²·n(θ), Kani/Robert gives
  ψ_θ as an explicit isogeny, and ker φ̂ = ker ψ_θ ∩ E_C[d] (generic case,
  CV §3.2). So the torsion requirement is **N > d·√n(θ)**, and the receiver's
  cost is one 2^e-isogeny chain in dimension 2, 4 or 8.

Everything else in the hint taxonomy is either decodable by everyone
(unmasked, single point with square N), or not decodable by anyone known
(Borel; any hint with N² < d — the surveys found no End-only decoder there,
and no impossibility proof either).

## 4. Proposal: Eigen-Séta (working name)

Séta's trapdoor one-way function with the torsion images masked by a
sender-random **diagonal** matrix in a **secret eigenbasis**. Equivalently:
the Castryck–Vercauteren FESTA backdoor, used by its rightful owner.

### 4.1 Scheme (KEM form; derived, not yet implemented)

Parameters: p = 2^e·3^f·c − 1, N = 2^e, d = 3^f, with 3^f ≈ 2^{2λ} and e
fixed by §4.3. Everything lives in F_{p²}.

**KeyGen.**
1. Walk from E_0 (j = 1728) with a secret isogeny of degree > 2^{4λ} to a
   statistically uniform E_A; keep the connecting ideal, so End(E_A) = O_A is
   known explicitly (SQIsign-style; Deuring tooling: "Deuring for the people"
   2023/106, ideal-to-isogeny in dimension 2 as in SQIsign2D-West 2024/760).
2. Enumerate short non-scalar θ ∈ O_A (LLL on the rank-3 trace-zero lattice;
   shortest norm ≈ p^{2/3} for a random maximal order, CV Remark 9). Pick θ
   with trace t, norm w, Δ' := 4w − t² such that
   (i) −Δ' is a non-zero square modulo 2^e (so θ is diagonalisable on E_A[2^e]
       with distinct eigenvalues μ ≠ ν mod 2), and
   (ii) 4N² − d²Δ' is a prime ≡ 1 mod 4 (for a dimension-2 decryption, §4.2).
   Re-sample θ (thousands of short vectors are available) or E_A until both
   hold; (i) has density ≈ 1/4 among odd Δ', (ii) ≈ 1/(2 ln(4N²)) ≈ 1/700.
3. Compute an eigenbasis (P,Q) of θ on E_A[2^e] (evaluate θ on a basis via its
   explicit representation; solve the 2×2 eigenproblem mod 2^e).
4. pk = (E_A, P, Q). sk = (θ as an evaluable endomorphism, μ, ν, t, w, and the
   Cornacchia solution of (ii)).

**Encaps.**
1. Sample s ← Z/3^f and nonces a, b ← (Z/2^e)^×, all derived from a seed.
2. K := P_d + [s]Q_d for a canonical basis (P_d, Q_d) of E_A[3^f]; φ: E_A →
   E_C := E_A/⟨K⟩ (3-isogeny chain, length f).
3. ct = (E_C, S := [a]φ(P), T := [b]φ(Q)); key = H(seed, ct).

**Decaps.**
1. Define ψ := φθφ̂ ∈ End(E_C), degree d²w. Then ψ(S) = [dμ]S and
   ψ(T) = [dν]T, and ψ̂ = φ(t − θ)φ̂ acts as [dν] on S and [dμ] on T.
   (Derived: ψ(aφP) = aφθ(φ̂φP) = a d φ(θP) = a d μ φ(P). The nonces never
   enter; no commutation argument is needed beyond this.) So ψ's action on
   E_C[2^e] is known for free.
2. Kani in dimension 2: with integers c, c' from the Cornacchia solution of
   (2c + dt)² + (2c')² = 4N² − d²Δ' (equivalently deg(ψ + [c]) + c'² = N²),
   the endomorphism F = [[ψ + c, c'], [−c', ψ̂ + c]] of E_C² is an N²-isogeny;
   compute it as two N-isogeny halves from the actions of ψ and ψ̂ on E_C[N]
   (standard SQIsignHD/FESTA technique), giving ψ as an evaluable map.
3. ker φ̂ = ker ψ ∩ E_C[3^f] (ψ is cyclic on the 3-part for generic θ);
   recover φ, hence K and s; recompute a, b from the seed, re-encrypt, compare
   (FO with implicit rejection; the eigenvalues μ, ν are secret and used in
   step 1, so malformed ciphertexts must never reach a non-rejecting path —
   compare the Moriya–Onuki adaptive attacks on FESTA).

### 4.2 Why dimension 2 is available without Séta-style curve engineering

Kani needs deg(ψ + [c]) + c'² = N². Expanding, this is
4N² − d²Δ' = (2c + dt)² + (2c')². The left side depends on θ only through Δ',
so the receiver *chooses* θ among its short endomorphisms until the left side
is a prime ≡ 1 mod 4 and solves with Cornacchia (parity of 2c + dt matches
dt automatically when the prime is odd). Séta had to construct the curve
around a norm equation; here the curve is random and only the choice of θ is
constrained. If (ii) is dropped, dimension 4 works with any θ (write
4N² − d²Δ' − (2c + dt)² as a sum of two squares over a free c), dimension 8
always works.

### 4.3 Parameter estimate at λ = 128 (derived, not measured)

Constraints: MITM on the 3^f-isogeny ⇒ 3^f ≳ 2^{2λ}; decoding needs
N > d√w; w ≈ p^{2/3} for a random curve (CV Remark 9); p ≈ N·d.
Solving N ≈ d·(Nd)^{1/3} gives N ≈ d², p ≈ d³:

| quantity | value | note |
|---|---|---|
| d = 3^f | 3^162 ≈ 2^257 | sender's degree |
| N = 2^e | e ≈ 515 | eigenbasis torsion; N > d·p^{1/3} |
| p | ≈ 2^772 | vs FESTA 1292, POKÉ 431, SILBE 13013 bits |
| pk | ≈ 2·97 B (E_A) + 2·65 B (P,Q compressed relative to a canonical basis) ≈ 325 B | |
| ct | ≈ 193 B (E_C) + ≈ 130 B (S,T compressed) ≈ 325 B | |
| decaps | one (2^515, 2^515)-isogeny chain in dimension 2 plus a 3^162 kernel recovery | FESTA-class cost; C implementations of comparable chains run in tens of ms (SQIsign2D-West verification class) |
| key recovery | EndRing(E_A): classical p^{1/2} ≈ 2^386, heuristic p^{1/3} ≈ 2^257, quantum p^{1/4} ≈ 2^193 | far above target; the ciphertext problem is the binding constraint |

If one instead engineers θ Séta-style with w ≈ 4d² (so that 4N² − d²Δ' is
small), the same p ≈ d³ results; the random-curve route is preferred because
its key-recovery problem is plain EndRing plus a leakage term (§5.1).

### 4.4 What is new relative to each neighbour

* vs **Séta**: masks (diagonal, sender-random) make the hint useless to
  Robert's generic attack; torsion requirement N > d√w instead of N > D²;
  no trapdoor-curve construction, random curve plus a short endomorphism.
* vs **SILBE**: diagonal instead of scalar masks removes the ω(N) ≥ 2λ
  constraint (13 013 → ≈ 772 bits), no F_p-rational E_0, no Frobenius; a
  shortest endomorphism of E_A replaces the degree-√p connecting isogeny,
  which also drops the lollipop degree from (d_A d_B)² to d²·p^{2/3}.
  Lost: SILBE's pk is *any* curve with no extra data; ours leaks an
  eigenbasis (§5.1).
* vs **FESTA/QFESTA/POKÉ**: pk carries no torsion images of a receiver
  isogeny from a known-End curve, so key recovery is an EndRing-type problem
  rather than CIST/C-POKE; no E_0 at encryption time at all.
* vs **CV §5.3**: the "backdoor" is now the honest receiver's secret, so the
  basis-hashing countermeasure does not apply.

## 5. Cryptanalysis checklist (what must be settled before anything is claimed)

### 5.1 Key recovery: "EndRing with an eigenbasis hint"
Problem: given (E_A, P, Q) with (P,Q) an eigenbasis of some θ ∈ End(E_A) of
norm w ≈ p^{2/3}, find θ (then End(E_A) follows at cost |disc|^{1/4} ≈ p^{1/6}
classically / subexponential quantumly via HLW 2023/1448).
Known approaches and their cost:
* Guess the eigenvalues modulo L' | 2^e with L'² > w, then Kani: ≈ L'² ≈ w ≈
  2^{4λ} guesses; no per-prime localisation because N is a 2-power and the
  (2,2)-chain fails only at the end (to be verified experimentally, §6).
* Cycle/MITM search for an endomorphism of norm ≈ w: √w ≈ 2^{2λ} classical,
  quantum claw w^{1/3} ≈ 2^{4λ/3}.
* The eigenlines ⟨P⟩, ⟨Q⟩ expose two horizontal 2-isogeny walks of length e
  in the Z[θ]-orientation (E_A/⟨2^k P⟩ are all Z[θ]-oriented). This is
  class-group information of length e ≪ h(Z[θ]) ≈ √w ≈ 2^{2λ}; the
  self-pairing attacks (Castryck et al. CRYPTO 2023, 2023/549) need m² | Δ'
  with m | N, excluded by choosing Δ' odd. Macula–Stange sesquilinear
  pairings (AC 2024) must be checked against the eigenbasis specifically.
* Distinguishing an eigenbasis of a norm-p^{2/3} endomorphism from a random
  basis: every basis is an eigenbasis of *some* endomorphism of norm
  O(p^{2/3} N^{4/3}) (CV §5.4), so the hidden quantity is only the *size*
  of the smallest diagonal endomorphism. State as a decisional assumption.

### 5.2 Message recovery: "eigen-CIST"
Given (E_A, P, Q, E_C, [a]φP, [b]φQ) with a, b random, find φ.
* Weil pairing leaks ab·d mod 2^e, so ≈ 2^{e−1} masks remain; the attacker
  needs correct masks on 2^{e'} with 2^{2e'} > d, cost ≈ 2^{e'} ≈ √d ≈ 2^λ.
  This is FESTA's CIST sizing (FESTA uses d of 257–260 bits at λ = 128).
* The exposed subgroups φ(⟨P⟩), φ(⟨Q⟩) (mask-independent) give SIDH-type
  4-tuples without torsion, believed hard; but these are *horizontal* in an
  unknown orientation and the attacker knows E_C/φ(⟨P⟩) ≅ pushforward of the
  oriented neighbour — check whether the Arpin–Castryck–Eriksen–Lorenzon–
  Vercauteren level-structure action (2024/1172) gives any lever.
* Level-structure reductions of De Feo–Fouotsa–Panny (2024/459, Cor. 12–13)
  must be re-run with the eigenbasis as the level structure.
* Leakage attacks on masks (Das–Invernizzi–Kutas–Meers PKC 2026, 2025/2075;
  fault attacks 2026/1944) apply as to FESTA; nonces must be fresh per
  ciphertext and derived from the seed.

### 5.3 Quantum
* Unstructured key recovery stays at p^{1/4}; the eigenbasis gives at best
  the §5.1 costs. Nothing here is a group action, so Kuperberg does not
  apply to the attacker (the orientation Z[θ] is hidden).
* If θ ever leaks (side channel on decaps), security drops to HLW's
  subexponential bound in |Δ'| ≈ p^{2/3}: this is the scheme's single point
  of failure and argues for constant-time Kani code and FO with implicit
  rejection from day one.

### 5.4 Correctness corner cases
* ψ non-cyclic on E_C[3^f] (ker φ θ-stable): probability ≈ 1/3 per prime
  factor, handled by CV §3.2's guessing (≤ 2 choices); with d = 3^f a
  single prime, expect O(1) retries or reject such K at encryption time via
  a public test (none exists without θ; so handle on the receiver side).
* Eigenvalues must satisfy μ ≢ ν mod 2 for the eigenlines to span E_A[2^e];
  guaranteed when Δ' is odd.

## 6. Experiment protocol (next steps, in order)

All runs on public toy parameters; no security claims at toy sizes.

1. **Correctness prototype (SageMath).** Toy p = 2^e·3^f·c − 1 with e ≈ 60,
   f ≈ 20. Use E_0 = j 1728 directly as E_A (End public is fine for
   correctness), θ = a short element of ⟨1, i, (i+π)/2, (1+iπ)/2⟩ with odd Δ'
   and −Δ' a square mod 2^e; eigenbasis; encaps with random diagonal masks;
   decaps via the dimension-2 (2^e, 2^e)-chain (reuse the FESTA-SageMath or
   Theta-Isogenies chain code) and the two-halves trick. Acceptance: decaps
   returns s for 1000 random seeds; measure the ψ-non-cyclic failure rate.
2. **Random-curve keygen.** Replace E_0 by a walked E_A with explicit O_A
   (SQIsign2D-West / "Deuring for the people" tooling); LLL for short θ;
   measure how many θ must be tried for conditions (i)+(ii) and the achieved
   w/p^{2/3} ratio. Acceptance: median w ≤ 4·p^{2/3}.
3. **Toy cryptanalysis of §5.1.** At e ≈ 40: brute-force eigenvalue guessing
   plus Kani, and record whether a wrong guess is detectable *before* the end
   of the chain (gluing/splitting pattern of intermediate varieties). If an
   early-abort signal exists the guessing cost drops from L'² toward e·4 and
   the design is dead; this is the most important single experiment.
4. **Toy cryptanalysis of §5.2.** Run the published Kani attack code against
   ciphertexts with (a) correct, (b) wrong masks; confirm failure mode is the
   final non-splitting only. Then the level-structure reductions.
5. **Parameter script.** Encode the constraints of §4.3 (MITM, N > d√w,
   p ≈ Nd, 2^e | p+1, 3^f | p+1, p^{1/3} and p^{1/4} margins) and emit
   (e, f, c) for λ ∈ {128, 192, 256}; compare bytes and chain lengths with
   FESTA, QFESTA, POKÉ, SILBE, FESTA-UPKE.
6. **Write-up target.** If 3 and 4 pass: define Eigen-EndRing and eigen-CIST
   formally, prove OW-CPA ⇐ eigen-CIST and key recovery ⇐ Eigen-EndRing,
   apply FO, and position against 2026/1431 Problem 7 (it would *not* solve
   Problem 7 — the ciphertext assumption is still a torsion-hint assumption —
   but it would be the first practical-size PKE whose key recovery is an
   EndRing variant with a leakage term rather than CIST/C-POKE).

## 7. Alternatives considered and why they are not the lead

* **Scalar-mask lollipop with a short secret endomorphism** (SILBE minus
  Frobenius): same ω(N) ≥ 2λ wall; p stays > 10 000 bits. Dropped.
* **Secret or partially published orientation + class group action** (large
  discriminant, PEGASIS-evaluable): the sender must act, so the action and
  hence Kuperberg are available to the attacker; OSIDH's fate. Dropped for
  the user's "EndRing-based" goal; kept as a benchmark.
* **Deuring Diffie–Hellman** (both parties know End of their own curve):
  no canonical function of (E_A, E_B) computable from either End alone is
  known; non-commutativity is exactly the obstruction. Research question,
  not a design.
* **Insufficient-torsion hints (N² < d) decodable with End**: no decoder and
  no impossibility result exists (surveys, §5 of the torsion-attack report).
  Genuinely open; a positive result here would be a new mechanism. Listed as
  the long-shot follow-up.
* **Suborder/IsERP-based trapdoors** (pSIDH-like): quantum-polynomial break.
  Excluded.
* **IBE from End(E)**: shown unrealisable or insecure (Ozbay Gurler–Struck
  SAC 2025, 2025/1726). Excluded.

## 8. Requirement status

| Requirement (from the request) | Status | Evidence / gap |
|---|---|---|
| Research the state of the art on EndRing-based encryption | complete (survey level) | Sections 1–2, sources below; ePrint PDFs were Cloudflare-blocked, mirrors used; items marked abstract-only in the agent reports were not re-verified from full text |
| Identify primitives/concepts to build on | complete | §3 decoding lemma (CV Lemma 3), Deuring tooling, dimension-2 Kani, FO |
| Devise a strategy for a novel scheme | complete as a design; **not implemented, not analysed beyond §5** | §4–6; novelty claim rests on the four surveys finding no eigenbasis-trapdoor PKE; a dedicated search for the exact combination is still advisable before any write-up |
| Parameter/efficiency figures | derived estimates only | §4.3; no measurements exist |
| Security | **not established** | §5 checklist; experiment 3 is the kill-or-keep test |

## Sources (primary; see agent reports for the full list)

* Wesolowski, FOCS 2021, ePrint 2021/919. Page–Wesolowski, EC 2024, 2023/1399.
  Herlédan Le Merdy–Wesolowski, TCC 2025, 2025/271; CiC 2025, 2023/1448.
* Wesolowski, 2026/1486; Mamah, 2026/1821; Udovenko, 2026/1575.
* Castryck–De Feo–Galbraith–Kutas–Reijnders–Wesolowski, "The isogeny problems", 2026/1431.
* Castryck–Vercauteren, AC 2023, 2023/1433 (Lirias full text).
* Fouotsa–Moriya–Petit, EC 2023, 2023/013. de Quehen et al., Crypto 2021, 2020/633.
  Petit, AC 2017, 2017/571. Sica, Math. Cryptology 2024.
* De Feo et al., Séta, AC 2021, 2019/1291 (HAL full text).
* Duparc–Fouotsa–Vaudenay, SILBE, SAC 2024, 2024/400 (EPFL full text, SAC slides).
* Basso–Fouotsa–Kouider–Kutas–Maino–Marco, UPKE from FESTA, AC 2026, 2026/1014.
* Basso–Maino–Pope, FESTA, AC 2023, 2023/660. Nakagawa–Onuki, QFESTA, Crypto 2024,
  2023/1468. Basso–Maino, POKÉ, EC 2025, 2024/624. PIKE 2026/473; INKE 2025/1458.
* De Feo–Fouotsa–Panny, level structure, EC 2024, 2024/459.
  Arpin–Castryck–Eriksen–Lorenzon–Vercauteren, 2024/1172.
* Chen–Imran–Ivanyos–Kutas–Leroux–Petit, AC 2023, 2023/779; Chen–Petit, EC 2025.
* Castryck–Houben–Merz–Mula–van Buuren–Vercauteren, Crypto 2023, 2023/549.
* Robert, EC 2023, 2022/1038. Page–Robert, Clapoti, 2023/1766. PEGASIS 2025/401;
  qt-Pegasis 2025/1859; CORAL 2026/896; Bonnetain et al. 2026/2215.
* Dartois–De Feo, PKC 2022, 2021/1681. Ozbay Gurler–Struck, SAC 2025, 2025/1726.
* Das–Invernizzi–Kutas–Meers, PKC 2026, 2025/2075. Gilchrist–Lai–Meyer, AC 2026, 2026/1944.
* SQIsign v2 specification (2025-07-07), §10.2.
