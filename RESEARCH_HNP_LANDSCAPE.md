# Side-Channel / HNP Attack Landscape for ECC

A focused survey of the **scalar-side** attack family on ECC — Hidden-
Number Problem (HNP) variants, lattice attacks on biased nonces,
Bleichenbacher signature attacks — and where it sits relative to
the **group-side** isogeny-graph attacks programme documented in
[`RESEARCH_SECP256K1_CM.md`](RESEARCH_SECP256K1_CM.md) and
[`RESEARCH_VOLCANO_FLOOR_RHO.md`](RESEARCH_VOLCANO_FLOOR_RHO.md).

**TL;DR**: side-channel/HNP attacks attack the *protocol layer*
(ECDSA, ECDH), not the *group structure*.  They are orthogonal to
isogeny-graph cryptanalysis.  They are also where every recent
real-world break of secp256k1/ECDSA has lived — TLS implementations
leaking nonce bits, hardware tokens with biased RNGs, fault attacks
on signing.  The isogeny-graph programme produced a structural-
completeness null result; the HNP programme produces *real attacks*
on real implementations.

## 1. The two attack surfaces

| Surface         | What is attacked                            | Hard problem                        | Mechanism                            |
|-----------------|---------------------------------------------|--------------------------------------|--------------------------------------|
| **Group-side**  | The mathematical group `E(F_p)`             | ECDLP                                | ρ, isogeny graphs, pairings, covers  |
| **Scalar-side** | Knowledge of the *secret scalar* `d` itself | Find `d` from partial info           | HNP/lattice, side-channel, fault     |

The group-side programme assumes a *cleanly-implemented*
black-box `(P, Q = d·P)` oracle and asks: can the group's
algebraic structure leak `d`?  The structural-completeness statement
we proved: **no, for prime-field CM curves, ρ is optimal up to
constant factors**.

The scalar-side programme assumes a *leaky implementation* and asks:
can partial information about `d` (or `k` in signatures) be
amplified into full recovery via lattice methods?  The answer is
**yes, almost always, given enough leakage**.

These are different attacks against different threats; not
substitutable.

## 2. The Hidden Number Problem (HNP)

**Definition** (Boneh–Venkatesan 1996): given a fixed unknown `α ∈
Z/nZ` and a collection of pairs `(t_i, a_i)` with

```
        a_i  =  ⌊α · t_i mod n⌋_{\text{high bits}}  +  e_i      (i = 1, …, m)
```

where `e_i` is the "noise" (unknown low bits), recover `α`.

If the high bits leak `k` bits per query and we have `m` queries,
HNP is solvable via lattice basis reduction (LLL/BKZ) when

```
        m · k  ≳  log₂(n)  +  small const.
```

i.e., total leaked bits exceed the secret size.  Typical:
3 bits/query × 100 queries for a 256-bit secret.

## 3. HNP-in-ECDSA: the canonical instantiation

ECDSA signing produces `(r, s)` from a random nonce `k`:

```
        r = (k · G).x  mod n
        s = k⁻¹ · (h + d · r)  mod n
        ⇒  k = s⁻¹ · (h + d · r)  mod n
```

If `k` is BIASED — even a few high bits known (or zero, or fixed) —
then for each signature `i`:

```
        d · r_i ≡  s_i · k_i  −  h_i   (mod n)
        ⇒  d · (r_i · s_i⁻¹) ≡  k_i  −  (h_i · s_i⁻¹)   (mod n)
```

This has the **exact HNP shape** with `α = d`, `t_i = r_i · s_i⁻¹`,
`a_i = bias-of-k_i`.  Lattice reduction on `m ≈ log₂(n) / bias-bits`
signatures recovers `d`.

### Real-world failures via this exact mechanism

- **Sony PS3 (2010)**: same `k` reused for every signature.  Trivial
  recovery: `d = (s_1 - s_2)⁻¹ · (h_1 - h_2) - r⁻¹ · s_1·(...)`.
  No lattice needed.
- **Android Bitcoin wallet (2013)**: bad `SecureRandom`, biased `k`.
  Multiple Bitcoin private keys recovered by passive observation.
- **TPM-FAIL (2019)**: timing leak of `k` bits during ECDSA signing
  in Intel/STMicro TPMs.  Allowed key recovery in minutes.
- **Minerva (2019)**: timing leak via scalar-multiplication length.
- **LadderLeak (2020)**: lattice attack on the Montgomery ladder
  used by some libraries, recovering keys from ~1 bit/query × `n`
  queries.

In all cases the attack worked on the **scalar side**, exploiting
implementation flaws.  The group structure of secp256k1 (j=0, CM,
sextic twists) was **irrelevant** to the attacks.

## 4. The codebase's HNP modules

| Module                                       | What it implements                                                                  |
|----------------------------------------------|-------------------------------------------------------------------------------------|
| `crypto/src/cryptanalysis/hnp_ecdsa.rs`      | Standard HNP attack on biased ECDSA nonces (Boneh–Venkatesan + Nguyen–Shparlinski). |
| `crypto/src/cryptanalysis/multi_key_hnp.rs`  | Multi-key HNP (BBSP): a single biased-nonce stream signed for multiple keys, allowing simultaneous recovery via a stacked lattice. |
| `crypto/src/cryptanalysis/bleichenbacher.rs` | Bleichenbacher's biased-nonce attack via DFT — works when bias is "fractional bits" too small for LLL. |
| `crypto/src/cryptanalysis/ecdsa_audit.rs`    | Signature-corpus auditor: detects repeated nonces, biased-nonce statistical signatures, fault patterns. |

These are the most practically-attack-relevant modules in the
codebase.  Every working ECDSA-breakable real-world incident in the
past 15 years is in their threat model.

## 5. Where HNP and isogeny-graph attacks intersect (they don't)

A common framing question: could the rich group structure of
secp256k1 (CM, sextic twists, isogeny graph) **amplify** an HNP
attack — i.e., make fewer biased signatures suffice?

The answer is **no, in any structural sense**, for several reasons:

1. **HNP works on `Z/nZ`**, not on `E(F_p)`.  The structure of `n` as
   a prime integer is what matters; the group structure of `E(F_p)`
   is irrelevant.

2. **Isogeny transport preserves `d`**.  If we transport
   `(P, Q = d·P)` to an isogenous `E'`, the same `d` appears on the
   new curve.  HNP signatures on `E'` give the same `d` mod `n`.
   No new information.

3. **GLV decomposition doesn't reduce HNP cost**.  Writing
   `d = d_1 + d_2 · λ` with smaller `d_1, d_2` does not help HNP,
   because the LATTICE problem still needs to recover `d` modulo
   `n`, not just `d_1, d_2` modulo `n/?`.  GLV is a multiplication
   accelerator, not a key-recovery accelerator.

4. **No published HNP variant uses j-invariant or isogeny data**.
   I checked: every HNP-on-ECDSA paper I'm aware of (Boneh–
   Venkatesan, Howgrave-Graham–Smart, Nguyen–Shparlinski,
   Bleichenbacher, more recent: Mulder et al, De Mulder et al,
   Aranha–Castaneda–Rivera, BBSP, Minerva) works *purely* on `Z/nZ`
   and lattice arithmetic.  The curve appears only via `r_i` in the
   signature.

The two attack programmes are **independent threat models** with
no known interaction.

## 6. The actual scalar-side open questions

Where would HNP research push the frontier?

- **Sub-bit-leak attacks**.  Bleichenbacher attacks via DFT
  recover `d` from sub-bit (fractional) bias in `k`.  Current
  state: ≈ 0.1 bits/sig with ~10⁹ signatures.  Pushing this below
  0.01 bits/sig would threaten more implementations.
- **Implicit-rejection HNP**.  When ECDSA signatures occasionally
  fail (e.g., `r = 0`), the failure pattern leaks `k` info.
  Quantifying this is open.
- **Multi-target HNP for many keys**.  BBSP-style lattices for
  10⁶+ keys simultaneously — feasibility at scale unclear.
- **Post-quantum-aware HNP**.  Grover/QFT speedups on the lattice
  reduction step — quantum HNP feasibility.

These are real research questions with real attack implications.
They are *out of scope* for the isogeny-graph cryptanalysis
programme.

## 7. Implementation defenses that make HNP useless

The HNP-cure stack:

1. **Deterministic ECDSA** (RFC 6979): `k = HMAC-DRBG(d, h)` — no
   randomness, no nonce bias.  Adopted by Bitcoin (libsecp256k1),
   most modern signing implementations.
2. **Ed25519** (Bernstein 2011): `k = SHA-512(d || h)[: half]` —
   structurally deterministic.
3. **Hedged nonces** (Aranha–Novaes–Takahashi 2020): `k =
   PRF(d, h, fresh-randomness)` — robust to both bias *and* RNG
   compromise.
4. **Implementation auditing**: corpora of ECDSA signatures
   statistically tested for nonce bias.  `crypto/src/cryptanalysis/
   ecdsa_audit.rs` is one such tool.

These defenses **completely close** the HNP attack surface for
correctly-deployed signing.  Where attacks still happen: bad RNGs
(JavaScript on a vending machine), hardware tokens that ignore RFC
6979, fault injection that violates the determinism.

## 8. Where the isogeny-graph and HNP programmes touch

There is **one** narrow scenario where the two interact:

> If an attacker has both (a) a biased-nonce HNP leak on ECDSA and
> (b) the ability to force signing on a *twist* of secp256k1 via an
> invalid-curve attack, the biased-nonce information leaks
> information about `d` *projected* into the twist subgroup.

This is a strict subset of the existing invalid-curve attack model,
not a new attack family.  The twist-subgroup leak is bounded by the
twist-order smooth factor (≤ 18 bits for secp256k1's worst twist,
per `RESEARCH_SECP256K1_CM.md` §2), so even with HNP combined with
invalid-curve access, the total leak is bounded by the twist's
smooth-part — no new "isogeny-graph enables HNP" amplification.

## 9. Verdict

The HNP/side-channel attack programme is **where real ECC breaks
happen**.  It is a **completely different attack surface** from the
isogeny-graph programme.  The structural-completeness statement we
proved (no isogeny-graph attack beats ρ for prime-field ECC) is
**unaffected by HNP advances**.  Conversely, advances in HNP do
not affect the isogeny-graph structural-completeness statement.

For an organisation deploying secp256k1, the threat ordering is
roughly:

1. **Bad RNG** in nonce generation (huge, real)
2. **Side-channel leakage** of `k` bits (real, attackable per case)
3. **Implementation bugs** (signature malleability, parsing, etc.)
4. **Key management failures** (private keys leaked via process)
5. ...
6. **ECDLP advances** (theoretical, no real progress in 30+ years)
7. **Isogeny-graph cryptanalysis** (structural-completeness null,
   per the programme in this repo)

For an organisation researching attacks, the value gradient is
roughly the reverse: structural-completeness null results are the
*easiest* to publish (well-defined, finite to verify), while real-
world implementation attacks are case-by-case empirical work.

## 10. References

- Boneh, Venkatesan, "Hardness of computing the most significant
  bits of secret keys in Diffie–Hellman and related schemes",
  Crypto 1996.
- Howgrave-Graham, Smart, "Lattice attacks on digital signature
  schemes", DCC 2001.
- Bleichenbacher 2002, talk: "A non-uniformity-based attack on
  DSA-SHA1".
- Nguyen, Shparlinski, "The Insecurity of the Digital Signature
  Algorithm with Partially Known Nonces", J. Cryptol. 2002.
- Bos, Halderman, Heninger, Moore, Naehrig, Wustrow, "Elliptic
  curve cryptography in practice", FC 2014. (Empirical survey of
  real ECDSA RNG failures.)
- TPM-FAIL paper, Moghimi et al, USENIX Security 2020.
- LadderLeak paper, Aldaya et al, CCS 2020.
- Aranha, Novaes, Takahashi, "Hedged ECDSA", J. Cryptol. 2020.
