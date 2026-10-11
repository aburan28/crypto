---

## 4.5 Best candidate algorithm: quotient-isogeny symmetrised index calculus (H1, with H4)

Scope: E/F_{qⁿ} (q = p or 2), rational torsion subgroup K = M^π of order k
(k ∈ {2, 3, 4}), target subgroup G of prime order r, h = N/r, decomposition
length m, subspace V ⊂ F_{qⁿ} of F_q-dimension ℓ_V. In the Koblitz case E/F₂
and K = ⟨(0, √b)⟩ are F₂-rational, so τ commutes with the quotient map.

```
Input : E/F_{q^n}, K = <T> rational torsion of order k, G = <P0> of order r,
        target Q0 in G, parameters m, V.
Output: log_{P0} Q0.

1. Quotient.      E2 <- E/K by Velu;  X <- x o phi  (rational function, degree k).
                  [cost: O(k) field ops once; X is stored as (num, den) polys]
2. Factor base.   F2 <- { Q in E2(F_{q^n}) : x(Q) in V }        (size ~ q^{l_V}/2)
                  F1 <- phi^{-1}(F2) ∩ E(F_{q^n})                (size ~ k |F2|)
                  Store for each Q in F2 ONE rational preimage P(Q) in F1
                  and the class table  F1 -> F2  (hash on x(phi(P))).
                  Unknowns: the logs of  [h]P(Q), one per Q in F2   (|F2| unknowns),
                  plus, Koblitz case, one per tau-orbit of F2        (|F2|/n).
                  [cost: |F2| evaluations of phi; memory |F2| (x, pointer)]
3. Relations.     repeat until  #relations >= |F2| (/n) + 20:
      3a. R <- a P0 + b Q0  (random a, b); Rq <- phi(R).
      3b. Solve  S_m^{E2}(u_1..u_m) with u_i in V   and the Weil-descent
          system in the symmetrised variables (elementary symmetric
          functions of u_1..u_{m}; the k-torsion symmetry is already
          absorbed because u = x o phi), by F4 / SAT.
          [cost c_D: dominated by F4 at the first fall degree; identical
           Macaulay shape to the plain S_m^{E2} system (Prop. 3.2)]
      3c. For each solution (u_i): Q_i <- the point of F2 with x = u_i
          (two sign choices each, resolved by one addition chain on E2);
          check  sum Q_i = Rq  on E2.
      3d. Lift:  P_i <- P(Q_i)  (table lookup, no search; Prop. 4.1).
          Emit the projected relation
                 a*log P0 + b*log Q0  =  sum_i log [h]P(Q_i)      (mod r)
          after multiplying both sides by h (precomputed  [h]P(Q) logs are
          the unknowns, so nothing further is needed).
          [cost: m table lookups + 1 verification on E2]
4. Linear algebra.  Sparse system over Z/r with |F2| (/n) columns, weight
                    <= m per row; Lanczos/Wiedemann.
5. Descent.         Express Q0 (or a second random R) via one more
                    decomposition; solve for log Q0.
```

**Data structures.** (i) A hash map x(φ(P)) ↦ (index of Q in F₂, one preimage
P). (ii) The symmetrised Macaulay matrix in the elementary symmetric functions
of the u_i restricted to V (the FHJRV "symmetrized" system). (iii) For the
Koblitz case, canonical τ-orbit representatives of F₂ (the existing Koblitz
ledger code already does this for F₁; the only change is that orbits are now
taken on E₂ and the class table replaces the point list).

**Complexity relative to the plain generator on E with factor base
{x ∈ V} of the same size |F₁| ≈ k|F₂|.** Collection: identical c_D per
target, identical success probability per target (the symmetrised system has
the same solution count as the plain system on E₂; by Thm. 4.2 each solution
packs k^{m−1} decompositions on E, all with the same projected logs), but the
number of required relations drops from |F₁| to |F₂| = |F₁|/k. Linear
algebra: (|F₁|/k)² instead of |F₁|². Net: ×k on collection, ×k² on linear
algebra, ×1 on everything else; preprocessing adds |F₂| isogeny evaluations
(negligible). Memory halves to quarters. This is the *entire* upside, and it
is exactly the upside the 2014–2015 papers measured.

**What would make it fail in practice.** If the implementation's F4 cost is
dominated not by the number of targets but by the first-fall degree (as the
binary ledger shows: the degree-4 Macaulay matrix is 96 % useful rows), then
collection cost per relation is unchanged and the only saving is the
×k fewer relations needed. The 2× milestone therefore requires k ≥ 2 and a
linear-algebra share that is not already negligible.
