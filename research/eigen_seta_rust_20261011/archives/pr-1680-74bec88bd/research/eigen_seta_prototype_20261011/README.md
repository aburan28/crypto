# Eigen-Séta toy prototype (2026-10-11)

SageMath prototype of the complete KeyGen / Encaps / Decaps path of the
"Eigen-Séta" design from `RESEARCH_ENDRING_PKE_STRATEGY.md` (section 4):
public key = curve plus an eigenbasis of a secret endomorphism on the 2-power
torsion; ciphertext = codomain of a 3^f-isogeny plus torsion images masked by a
sender-random diagonal matrix; decryption = lollipop endomorphism (masks cancel
on the eigenbasis) decoded with a (2^e, 2^e)-isogeny chain in dimension 2.

> **Correction (2026-10-11).** The `walk` variant below is **insecure**: its
> decryption key is computed from Δ' = 7·25^k alone (`eigen_seta_toy.py`
> lines 133–155) and k is public, so anyone can decrypt (an exploratory run,
> not archived, decrypted 5/5 ciphertexts this way). The "non-cyclic ψ'"
> corner case reported below was never explained; it cannot occur in the
> secret-degree v2 scheme (`papers/eigen-seta/`, Proposition 3). The PDF
> `eigen_seta_design_note.pdf` predates these corrections. The maintained
> implementation is the Rust crate `eigen_seta/` (v2, secret degree,
> constant-time decapsulation, NIST-level parameter sets).

**Status.** Correctness only. Toy parameters (32–66-bit primes), no security,
no timing claims. Two *special-key* variants are implemented (see below); the
random-curve variant of the note is **not** implemented because it needs the
two-halves (N² > deg) Kani evaluation, which requires dual Richelot chains.

## Files

| file | role |
|---|---|
| `eigen_seta_toy.py` | the scheme: parameter search, keygen (variants `cm`, `walk`), encaps, decaps, lemma check, driver |
| `kani_chain.py` | generic (2^e,2^e)-chain between elliptic products allowing product→product, gluing, Richelot and splitting steps in any order (new; needed because the eigen-structure forces product-type initial steps) |
| `richelot_aux.py` | vendored, unmodified, from Castryck-Decru-SageMath (MIT; `LICENSE.richelot_aux`): gluing, Richelot correspondence, splitting, divisor doubling |
| `wrong_guess_experiment.py` | protocol experiment 3: do wrong but pairing-consistent guesses of the secret discriminant fail early or only at the final split? |
| `results/` | raw logs of the runs quoted below |

## Running

Needs SageMath (tested with pip-installed `passagemath-standard` 10.8.12 on
Python 3.13; `python -c "from sage.all import *"` must work).

```
python eigen_seta_toy.py E F K TRIALS VARIANT
#   p = 2^(E+1) 3^F 5^K c - 1,  d = 3^F sender degree,  5^K secret walk degree
python eigen_seta_toy.py 26 5 0 25 cm      # CM(-7) curve as public key
python eigen_seta_toy.py 36 5 2 25 walk    # hidden CM curve behind a 5^2-walk
python wrong_guess_experiment.py 36 5 2 walk
```

Exit code 0 iff every decapsulation returned the sender's secret `s`.

## Variants

* `cm`: E_A is the curve with CM by Z[(1+√−7)/2] (j = −3375), θ = ω of degree 2
  (trace ±1, norm 2, Δ' = 7). End(E_A) is public, so this is insecure as a key;
  it exercises the whole mechanism with a tiny Kani chain.
* `walk`: E_A = ε(E_cm) for a secret 5^k-isogeny ε, θ = ε ω ε̂ (trace ±5^k,
  norm 2·25^k, Δ' = 7·25^k). This is the Castryck–Vercauteren §5.3 "backdoored
  FESTA basis" used by its legitimate owner. Chain condition 2^e > 7·25^k·9^f.

Both run the single-chain regime (torsion N > deg ψ'). The note's random-curve
variant needs N² > deg ψ' (two halves) and is the next milestone.

## What the runs showed (observations, toy sizes)

Raw logs: `results/correctness_runs.txt`, `results/wrong_guess_*.txt`.

* Decapsulation recovers `s` with the masks cancelling, in every completed
  trial, for both variants (cm: e = 21, 26; walk: e = 30, 36, 44; chain
  lengths 21–44; p from 32 to 66 bits). The independent `lemma_check`
  (direct evaluation of φθ'φ̂ on the ciphertext points) agrees with the
  predicted eigenvalues ±d·s₀ in every trial.
* The (2,2)-chain always starts with product-to-product steps (2–3 of them)
  before the gluing, and ends with product steps after the split, exactly as
  predicted by the eigen-structure (θ' ≡ 1 on E_C[2]).
* Corner cases seen: ker ψ' ∩ E_C[3^f] non-cyclic (handled by trying the
  cyclic order-d subgroups and keeping the one landing on E_A); a degenerate
  Mumford-coordinate input in the splitting formulas for some payload points
  (handled by retrying with a fresh basis of E_C[3^f]). Frequencies are in
  the logs; causes not isolated.
* Wrong-guess experiment (walk e=36: 1 true + 3 "low bits right" + 40 random
  consistent guesses; cm e=26: 1 + 3 + 40): every wrong guess ran the full
  chain and ended on a Jacobian (no split) — no early-abort signal — with one
  exception in the cm run where a random wrong Δ' also split to E_C × E_C. In
  the cm variant E_C has many small endomorphisms (it is 3^f-isogenous to a
  CM curve of discriminant −7), so other elements of the same quadratic order
  share the eigenbasis; this false positive is an artefact of the CM toy and
  is itself a warning for special-key variants.

## Not done

* Two-halves Kani evaluation (dual Richelot chain + middle-Jacobian
  isomorphism, or a theta-model implementation): required for random curves.
* Random-curve keygen (LLL on End(E_A) and evaluation of a short
  endomorphism): a first attempt via the kernel→ideal→right-order route showed
  that for curves 5^k-isogenous to E_0 every odd-trace endomorphism has norm
  ≥ p/(4·5^k), so that route cannot be exercised at toy size with one chain.
* Parameter script, timing, any security analysis beyond experiment 3.
