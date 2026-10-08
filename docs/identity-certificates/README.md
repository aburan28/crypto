# Identity certificates

**A claim that two polynomials are equal, as a record a stranger can
replay without trusting the author.**

Schema `identity.certificate/v1`, id `IDC1h<16 hex>`. Implemented natively
in `src/cryptanalysis/identity_certificate.rs`; the worked example lives in
this directory as `schoenhage-10mult.json` and is pinned by a unit test.

## 1. Why

Alman and Vassilevska Williams (`arXiv:2610.06783v1`) refuted the 3SUM and
APSP hypotheses with an algorithm whose base identity was found by machine
and whose main results were verified in Lean 4. This repository already
verifies **runs**: sealed `ecbench` sessions, replay receipts, bounds and
verdicts (`docs/bounds/`). It had no record for **identity** claims —
"polynomial `L` equals polynomial `R`" — which until now lived in prose, or
beside a unit test whose evaluation points the author chose. A checker that
trusts the author's points is one the author can satisfy by accident or by
construction. The next machine-found identity to reach this tree (a
bilinear identity, a summation-polynomial relation, an orbit-fold counting
argument) should enter through a record that nobody has to take on faith.

The study that proposed this format is
`research/lopsided_implications_20261006/` (method implication 1).

## 2. What a certificate is

| field | meaning |
|:--|:--|
| `schema` | always `identity.certificate/v1` |
| `id` | `IDC1h` + first 16 hex of SHA-256 over the **canonical statement bytes** (sorted-key compact JSON of `statement`) |
| `statement` | `name`, `source`, ordered `variables`, declared `degree_bound`, and two expression trees `lhs`, `rhs` (`var`, `const`, `neg`, `add`, `mul`) |
| `modulus` | `2^61 - 1` in v1 |
| `points` | number of evaluation points |
| `seed` | author-supplied; mixed with the `id` before any point is drawn |
| `value_digest` | SHA-256 over the little-endian left-hand values, in point order |
| `independence` | who issued and who checked, as text in the record |
| `proof` | optional `{kind, path, sha256}` pointer to a formal proof artifact |

Two tiers follow from the last field. **Evaluated** (`proof: null`): the
identity holds at every derived point over `F_p`, with the Schwartz–Zippel
bound below. **Proved** (`proof` set): additionally, a formal proof exists
at the named path with the named digest. The checker records that pointer;
it does not run the prover. A `proved` certificate is therefore a pointer
the reader must follow — exactly the status the repository's citation rules
give anything an agent has not itself opened.

## 3. What the checker does

`check(cert)`:

1. refuses a record whose `schema` or `modulus` is wrong, whose `points` is
   zero, whose variables are duplicated or undeclared, or whose declared
   `degree_bound` is below the **syntactic** degree of either side
   (`var = 1`, `const = 0`, `add = max`, `mul = sum`) — the bound is
   checked, not trusted;
2. recomputes the `id` from the statement and refuses on mismatch, so a
   statement cannot be swapped under a sealed id;
3. derives every point from SHA-256 of `(id ‖ seed ‖ "cert")` expanded by
   splitmix64, so the author cannot pick points independently of the
   statement being certified;
4. evaluates both sides at every point over `F_p` and **rejects on the
   first disagreement, returning the point as a counterexample**;
5. recomputes `value_digest` and rejects on mismatch;
6. with `check_with_extra(cert, (checker_seed, k))`, evaluates `k` further
   points from a domain-separated stream the certificate never saw.

`issue(statement, seed, points, independence, proof)` runs the same
evaluation and **refuses to seal a false identity**: a statement that fails
at any derived point returns the counterexample instead of a certificate.
A certificate for a false statement can therefore only exist by forging the
record, and forging the record is what step 4 catches.

**Soundness.** A nonzero polynomial of total degree at most `d` over `F_p`
vanishes at a uniformly random point with probability at most `d / p`, so
`k` independent points accept a false identity with probability at most
`(d / p)^k`; the verdict reports `log2` of that bound, `k · (log2 d − 61)`.
For the worked example (`d = 3`, `k = 32`) it is below `2^{−1900}`.

## 4. What a certificate does not prove

- It certifies the identity **over `F_p`, `p = 2^61 − 1`**. Two integer
  polynomials differing by a polynomial all of whose coefficients are
  divisible by `p` would pass. State coefficient magnitudes when the claim
  is meant over `Z`.
- It says nothing about the *usefulness* of the identity — rank, border
  rank, sparsity, or any algorithmic consequence. Those are separate claims
  with separate evidence; the Schoenhage certificate below says the
  identity holds, not that it yields `O(N^2 / D^0.063)`.
- The `proved` tier is a pointer. Nothing here runs Lean.

## 5. The worked example

`schoenhage-10mult.json` certifies Lemma 6 of the paper — Schoenhage's
ten-term identity computing a `3×3` outer product and a `2×2` inner product
together, plus a harmless error polynomial `E`, in 24 variables at degree 3.
It is the load-bearing identity under every speedup in that paper, which
makes it the right first thing to certify here. Regenerate with

```sh
cargo run --release --example identity_certificate_schoenhage \
  > docs/identity-certificates/schoenhage-10mult.json
```

The same seed reproduces the file byte for byte; any change to the
statement changes the id, which is the point. The unit test
`committed_schoenhage_certificate_still_checks` pins the committed file, and
`the_identity_without_its_error_term_is_false_and_cannot_be_issued` shows
the refusal path on the identity with `E` dropped — a false identity fails
at the first derived point and cannot be sealed.

## 6. When to write one

Any research round that claims an algebraic identity — in a report, a
module doc, or a comment that another result rests on — ships the
certificate beside the claim, in the study directory or here, and cites
its `IDC1h…` id. A claim without one is prose. This is the same rule the
repository applies to runs (a figure without its sealed session is prose)
extended to the one kind of claim it did not yet cover.
