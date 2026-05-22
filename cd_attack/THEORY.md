# The Castryck-Decru attack on SIDH — what it actually is

A self-contained explanation of the 2022 attack that broke SIDH/SIKE, anchored
to the implementation in this codebase. No prior knowledge of isogeny crypto
assumed; basic elliptic-curve fluency required.

---

## 1. The setting: SIDH and what it claimed

**SIDH** (Supersingular Isogeny Diffie-Hellman, De Feo–Jao–Plût 2011) is a
key-exchange protocol that ran on supersingular elliptic curves over
$\mathbb{F}_{p^2}$ for a prime of the form $p = 2^a \cdot 3^b - 1$. It looked
quantum-secure: it was an alternate-round NIST PQC finalist (SIKE) right up
until July 2022.

The setup:

- Two participants, **Alice** and **Bob**, share a public starting curve
  $E_0/\mathbb{F}_{p^2}$ and four public auxiliary points:
  - $P_A, Q_A \in E_0[2^a]$ — a basis of the $2^a$-torsion subgroup
  - $P_B, Q_B \in E_0[3^b]$ — a basis of the $3^b$-torsion subgroup

- **Alice's keygen**: pick secret $m_A \in [0, 2^a)$, form the kernel
  $R_A = P_A + m_A \cdot Q_A$, and walk a chain of $a$ degree-2 isogenies from
  $E_0$ killing this kernel, ending at her public curve $E_A$. She publishes:

  $$\boxed{E_A, \quad \varphi_A(P_B), \quad \varphi_A(Q_B)}$$

  where $\varphi_A: E_0 \to E_A$ is her secret isogeny (degree $2^a$).

- **Bob's keygen** is symmetric: pick $m_B \in [0, 3^b)$, kernel $R_B = P_B + m_B Q_B$, chain of $b$ degree-3 isogenies, publish $E_B$ and the images of Alice's basis.

- **Shared secret**: each party uses their secret to walk an isogeny *on the
  other party's published curve*, producing the same $j$-invariant.

The hardness assumption: even with the auxiliary points, recovering $m_A$ or
$m_B$ should take exponential time. The supersingular isogeny graph has
exponentially many vertices, and walking it without the kernel was believed
to be Pollard-rho-hard.

This was wrong.

## 2. The auxiliary points are far too generous

SIDH leaks more than people noticed. Knowing **just $E_A$**, recovering
$\varphi_A$ would indeed be a search through the isogeny graph. But the
published auxiliary points tell us exactly *where the isogeny sends two
specific points*:

$$\varphi_A(P_B), \quad \varphi_A(Q_B)$$

Since $\{P_B, Q_B\}$ is a basis of $E_0[3^b]$, this fully determines
$\varphi_A$'s action on the entire $3^b$-torsion. That's an enormous amount of
information for an isogeny of "size" $2^a$ — the isogeny is determined up to
isomorphism by its kernel (a single point of order $2^a$), and we get its
action on $2 \cdot 3^b$ different points.

The intuition for the attack is: this redundancy must be exploitable. The
question is *how*.

## 3. Kani's lemma (the load-bearing fact)

The key tool is a theorem from Ernst Kani, published 1997, originally about
Jacobians of curves. Modernized statement:

**Theorem (Kani's lemma).** Let $E$ be an elliptic curve, let $\varphi_1, \varphi_2: E \to E_1, E_2$ be isogenies of coprime degrees $d_1, d_2$ with $d_1 + d_2 = N$. Then the matrix

$$\Phi = \begin{pmatrix} \varphi_1 & \varphi_2 \\ -\widehat{\varphi_2} & \widehat{\varphi_1} \end{pmatrix} : E \times E \longrightarrow E_1 \times E_2$$

is an isogeny of abelian surfaces of degree $N^2$, with kernel structure
determined by how $\varphi_1, \varphi_2$ act on the $N$-torsion of $E$.

(Here $\widehat{\varphi}$ is the dual isogeny: it goes the other way, and
$\widehat{\varphi} \circ \varphi = [\deg \varphi]$.)

Why this is interesting: $\Phi$ is a $(N, N)$-isogeny — meaning its kernel
intersects each factor $E \times \{0\}$ and $\{0\} \times E$ in a subgroup of
order $N$. When $N$ is a prime power, **$\Phi$ can be decomposed as a chain
of $(\ell, \ell)$-isogenies on abelian surfaces**, computable explicitly.

The abelian surface $E \times E$ is generically a Jacobian of a genus-2 curve.
$(2, 2)$-isogenies on such Jacobians are called **Richelot isogenies** and
have explicit formulas (originally Richelot 1837, modernized by Bost-Mestre
1988). This is the mechanical engine.

## 4. The attack idea

Set $N = 2^a$ (Alice's torsion exponent). We want to find an isogeny
$\varphi_2: E_0 \to E_S$ of some degree $c$ with $c + 2^a = (\text{a power
of 2})$ — wait, that's not quite right. The actual relation is:

$$2^{\text{exp}} - 3^n = u^2 + 4 v^2$$

(see `src/uvtable.rs` and `tests::diophantine_relation_holds`)

The right-hand side $u^2 + 4v^2$ is the **norm** of the endomorphism
$\alpha = u + v \cdot (2\iota)$ on the special starting curve $E_0$, where
$\iota: E_{1728} \to E_{1728}$ is the j=1728 automorphism with $\iota^2 = [-1]$
(so $(2\iota)^2 = [-4]$, hence $\deg(2\iota) = 4$). On $E_0$, transported via
the 2-isogeny to $E_{1728}$, we get the same shape.

Now the construction: for each candidate guess $j$ of Bob's first $n$ ternary
digits, we form a *predicted* partial $3^n$-isogeny on $E_0$ with kernel

$$K_3 = 3^{b-n}(P_B + j \cdot Q_B)$$

If $j$ matches Bob's true digits, this $K_3$ generates a subgroup of Bob's
real kernel. We then form the **Kani-distorted kernel**:

$$K_\text{distort} = u \cdot K_3 + v \cdot (2\iota)(K_3)$$

This is the kernel of an isogeny $\tau$ on $E_0$ whose degree equals
$(u^2 + 4 v^2) \cdot 3^n = 2^{\text{exp}}$. Combined with what we know from
Bob's published data, **Kani's lemma's matrix $\Phi$ becomes computable** —
the row $(\varphi_1, \varphi_2) = (\text{Bob's } 3^n\text{-piece}, \tau)$
sums to a total degree that's a power of 2, exactly what we need for a
$(2, 2)$-chain on abelian surfaces.

Algorithmically:

1. Walk the predicted $\tau$ chain on $E_0$ to land on some new curve $C$.
2. Push Alice's $2^a$-torsion through this chain to get points $P_c, Q_c$ on $C$.
3. Use Bob's published auxiliary data $(\varphi_A(P_B), \varphi_A(Q_B))$ —
   actually wait, the directions are different: in the attack on Bob's
   secret, we use Alice's $2^a$ torsion images on Bob's curve $E_B$, namely
   $\varphi_B(P_A), \varphi_B(Q_A)$ — call these $P_B', Q_B'$.
4. **Glue** the elliptic product $C \times E_B$ into a single genus-2
   Jacobian $J(h)$ via Kani gluing.
5. Walk a chain of $\text{exp} - 2$ Richelot $(2,2)$-isogeny steps on $J(h)$,
   pushing $(P_c, P_B')$ and $(Q_c, Q_B')$ forward.
6. **Check the splitting**: does the chain end with a Jacobian that
   decomposes as a product of two elliptic curves?

**The recovery signal**: if our guess $j$ was correct, the Kani structure
is genuine and the chain *will split* at the end. If $j$ was wrong, the
chain either splits *prematurely* (somewhere in the middle, indicating
the construction is fictional) or doesn't split at all.

That's the whole attack: try each candidate $j$ in $\{0, 1, \ldots, 3^n - 1\}$, run the chain, see which $j$ gives a clean split. Recover the digits, advance, repeat.

The kicker: each "try" runs in polynomial time. The number of tries per
chunk is $3^n$ for chunk size $n$, but $n$ stays small (the uvtable gives
$n \le 193$ for SIKEp434, processed in chunks). Total recovery: hours on a
laptop for SIKEp434, where the security target was "centuries on a quantum
computer."

## 5. Why supersingular curves specifically

The construction needs a curve with **rich endomorphism structure** —
specifically, we need the integer $u + v(2\iota)$ to be an actual
endomorphism of $E_0$, not just an abstract algebraic combination. For an
*ordinary* elliptic curve, $\mathrm{End}(E)$ is an order in an imaginary
quadratic field, which has a small number of degree-bounded endomorphisms.

For a **supersingular** curve, $\mathrm{End}(E)$ is an order in a **quaternion
algebra** — a much richer object. The endomorphism $2\iota$ of $E_{1728}$ is
genuine and computable, and combined with $\mathbb{Z}$-multiples gives us
plenty of room to build Kani-style auxiliary isogenies.

SIDH chose supersingular curves precisely because their isogeny graph is
nicely "Ramanujan" (rapidly mixing), thinking that this would give better
security. The opposite is true: the rich endomorphism algebra is precisely
what makes the Kani construction work. This is a beautiful irony.

## 6. The Diophantine table

The uvtable (`src/uvtable.rs`) precomputes, for each odd chunk size $n$, the
smallest power exp such that

$$2^{\text{exp}} - 3^n = u^2 + 4 v^2$$

has a solution in non-negative integers $u, v$. The first entries:

| $n$ | $\exp$ | $u$ | $v$ | check |
|----:|------:|----:|----:|---|
| 1 | 3 | 1 | 1 | $8 - 3 = 5 = 1 + 4$ ✓ |
| 3 | 5 | 1 | 1 | $32 - 27 = 5 = 1 + 4$ ✓ |
| 5 | 8 | 3 | 1 | $256 - 243 = 13 = 9 + 4$ ✓ |
| 7 | 13 | 23 | 37 | $8192 - 2187 = 6005 = 529 + 5476$ ✓ |
| 9 | 15 | 59 | 49 | $32768 - 19683 = 13085 = 3481 + 9604$ ✓ |

The condition needs $2^{\text{exp}} - 3^n$ to be a sum of squares in a specific
way, which by Fermat's theorem on sums of two squares requires every prime
factor $\equiv 3 \pmod 4$ to appear to an even power. The uvtable is built
once per parameter family using number-theoretic search.

`tests::diophantine_relation_holds` in our code verifies this identity for
every table entry — runtime assertion that the table is internally consistent.

## 7. Concrete attack flow in this code

Mapping the theory to the implementation:

```
SIDH-side (Bob's public key, set up by `simulate_bob_keygen`):
  E_B = E_0 / ⟨P_B + m_B · Q_B⟩    [via pushing_3_chain, b steps]
  Published: (E_B, φ_B(P_A), φ_B(Q_A))

Attack-side (driven by `recover_bob_secret`):

for each chunk of n ternary digits to recover:
  for each candidate guess j ∈ [0, 3^n):
    # Build the "guessed Bob" 3^n-isogeny kernel
    K_3 = 3^(b−n) · (P_B + j · Q_B)              [on E_0]

    # Distort via Kani: add the auxiliary endomorphism u + v·(2ι)
    K_distort = u · K_3 + v · (2ι)(K_3)          [in `build_kani_source`]

    # Push the source through the predicted chain
    (C, P_c, Q_c) = Pushing3Chain(E_0, K_distort, n, [P_A, Q_A])

    # Glue C × E_B into a genus-2 Jacobian
    h = from_prod_to_jac(C.2-tor, E_B.2-tor)     [in `src/glue.rs`]

    # Compute the gluing image divisors via msolve Groebner basis
    D_1 = gluing_divisor_map(α, β, [2^alp·P_c, 2^alp·P_B'])   [via msolve_bridge.py]
    D_2 = gluing_divisor_map(α, β, [2^alp·Q_c, 2^alp·Q_B'])

    # Run the (2,2)-isogeny chain on J(h)
    split = does_22_chain_split(h, D_1, D_2, exp)   [in `src/richelot.rs`]

    if split:
      → j is the correct chunk of Bob's secret
      → accumulate, advance i by n, continue outer loop
```

The "splitting check" at each step is the determinant $\delta$ of the
partition matrix — see `src/richelot.rs::is_split`. Geometrically, $\delta = 0$
means the codomain of the next Richelot step is a product of two elliptic
curves (rather than a Jacobian of a smooth genus-2 curve), which is exactly
the Kani-end-state we expect for a correct guess.

## 8. Why this can't easily be patched

You can't simply "remove the auxiliary points" from SIDH — Bob *needs* them
to be able to compute the shared secret on his side (Alice's torsion images
on his curve $E_B$ are what let him pick the right kernel). The auxiliary
points are not a bug; they're the protocol's defining feature.

Variants:

- **SIKE** (the NIST submission) is identical structure and dies the same way.
- **CSIDH** (commutative isogeny DH) uses *only* the image curve, no
  auxiliary points. There's no Kani-lemma hook. CSIDH survives.
- **M-SIDH, MD-SIDH** (later "masked" variants) try to hide the auxiliary
  points behind quadratic-twist randomization. Several refinements have been
  proposed; none have full consensus on security.

## 9. Tying it back to the codebase

Every concept above has a concrete implementation:

| Concept | Code |
|---|---|
| Supersingular $E_0 = M_6$ over $\mathbb{F}_{p^2}$ | `ec::Montgomery` |
| The endomorphism $\iota, 2\iota$ | `ec::iota_on_e1728`, `ec::two_iota_on_e1728` |
| 3-isogeny chain ($\varphi_B$ or guessed kernel) | `ec::pushing_3_chain` |
| Kani-distorted source-side kernel | `attack::build_kani_source` |
| Kani gluing $E_\alpha \times E_\beta \to J(h)$ | `glue::from_prod_to_jac` |
| Genus-2 Jacobian arithmetic (Cantor) | `jacobian::add`, `jacobian::scalar_mul` |
| Richelot $(2,2)$-isogeny codomain | `richelot::codomain` |
| Richelot splitting test ($\delta = 0$) | `richelot::delta`, `richelot::is_split` |
| Richelot divisor pushforward | `richelot::pushforward` |
| Chain of Richelot steps | `richelot::does_22_chain_split` |
| The Diophantine Kani recipe | `uvtable::UVTABLE` |
| Per-candidate try loop | `attack::try_guesses` |
| Outer bit-recovery loop | `attack::recover_bob_secret` |
| Gluing-divisor-map via Groebner | `scripts/msolve_bridge.py` |

Each of the 42 unit tests checks a specific mathematical property — Mumford
validity after Cantor reduction, $|J(h)| = |E_\alpha| \cdot |E_\beta|$ for
the gluing, $\iota^2 = [-1]$, the Diophantine relation, etc. The tests are
not just "code runs" — each verifies a real algebraic invariant.

## 10. The historical surprise

Kani's lemma was published 1997. Bost-Mestre's explicit Richelot formulas
were 1988. The Castryck-Decru attack on SIDH was 2022 — a quarter century
after all the ingredients existed.

The puzzle isn't that the attack is hard once you know to look. It's that
nobody connected those tools to the SIDH security assumption. SIKE had been
through five years of intense scrutiny by professional cryptographers; the
attack took an algebraic geometer (Wouter Castryck) and a graduate student
(Thomas Decru) two months to find, after they decided to look.

The lesson: cryptographic security claims about novel structures benefit
heavily from cross-pollination with adjacent areas of pure math. The
*existence of auxiliary structure that the security proof ignores* is a
recurring source of cryptanalytic surprise — see also Frey-Rück on
hyperelliptic-curve DLP, Joux on pairing-based identity, the BGN
breakage of Boneh's identity-based encryption, etc.

---

## Recommended reading

- **Castryck, Decru. "An efficient key recovery attack on SIDH"**
  *Eurocrypt 2023* / IACR ePrint 2022/975. The original.
- **Maino, Martindale et al. "A direct key recovery attack on SIDH"**
  *Eurocrypt 2023* / IACR ePrint 2022/1026. Independent parallel discovery.
- **Robert. "Breaking SIDH in polynomial time"**
  *Eurocrypt 2023* / IACR ePrint 2022/1038. Generalization to non-CRT-friendly cases.
- **Kani. "The number of curves of genus two with elliptic differentials"**
  *J. Reine Angew. Math.* 1997. The original Kani's lemma paper.
- **Costello. "B-SIDH"** (and his SIDH tutorial). Background reading on
  what SIDH was supposed to be.

For the implementation side, the de facto reference is the
[Castryck-Decru Magma artifact](https://homes.esat.kuleuven.be/~wcastryc/),
which is ~250 lines of Magma code leaning heavily on Magma's built-in
genus-2 Jacobian and Groebner basis routines. The codebase here is a
from-scratch Rust reconstruction of those primitives, with msolve in place
of Magma for the Groebner step. Same math, much more line-auditable.
