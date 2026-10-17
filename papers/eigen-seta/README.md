# Eigen-Séta: an isogeny KEM whose trapdoor is a hidden-degree endomorphism

`eigen_seta.tex` specifies Eigen-Séta v2, proves that it is correct, derives
parameters for λ = 128, 192 and 256, analyses the attacks we know, and
describes and measures the Rust implementation in `crypto:eigen_seta/`.

**Status.** Research draft, anonymized, not peer reviewed, not submitted and
not posted. Security is a conjecture supported only by the analysis in the
paper. Restore the authors with `make deanonymize` after creating an
`AUTHORS` file (see `../AUTHORS.example`).

## What it says

* **Specification.** Key generation, encapsulation, decapsulation and the FO
  transform, exactly as implemented (§3).
* **Correctness.** Proofs of correctness (§4). Proposition 3 shows that the
  kernel recovered during decryption is always cyclic, so the corner case
  reported for the toy prototype cannot occur in v2.
* **Security.** Key recovery reduces to finding the secret degree f. The
  binding attack is message recovery by guessing the mask ratio modulo about
  √(3^F) and then running Robert's attack, at about 2^λ classical and
  2^(λ/2) Grover iterations. On that basis we *conjecture* NIST categories
  1, 3 and 5 (§6 and its attack-cost table).
* **Breaks of earlier variants.** The public-degree "walk" variant of the
  earlier toy prototype is insecure (§6.2).
* **Implementation.** Constant-time decapsulation, checked by instruction
  counts and a dudect-style test. The probe found a compiler-introduced
  branch, which is fixed (§7).
* **Measurements.** Sizes, timings, correctness rates and constant-time
  evidence (§8). Each number cites a committed record under
  `crypto:research/eigen_seta_rust_20261011/`.

## What it does not say

* No proof of security, and no claim that the scheme resists attacks we did
  not consider. The key assumption, that the eigenbasis gives no handle on f,
  is untested by anyone else.
* No performance claim beyond the measured host. No KEM-level speedup is
  claimed, because no baseline binary was kept.
* Key generation is not constant time. Memory-address constant time is
  argued, not measured. Power, EM and fault attacks are out of scope.

## Building

```
make           # PDF (pdflatex, bibtex, two more passes)
make check     # unresolved refs/citations, LaTeX errors, overfull boxes, page count
make clean
```

Figures are drawn from committed run records only, with:

```
python figures/plot_paper.py ../../research/eigen_seta_rust_20261011
```

## Reproduction

| What | Where |
|---|---|
| code, tests, KATs | `crypto:eigen_seta/` (`cargo test --release`) |
| parameters | `cargo run --release --example gen_params` |
| timings | `crypto:research/eigen_seta_rust_20261011/reports/run-20261011T1813Z-bench/` |
| correctness rates, constant-time probe at Toy16 and L1, dudect, f-filter rates | `crypto:research/eigen_seta_rust_20261011/reports/run-20261011T1955Z-verify/` |
| earlier constant-time probe (commit 19c641d12) | `crypto:research/eigen_seta_rust_20261011/reports/run-20261011T1827Z-sidechannel/` |

## Before submission

1. **Novelty audit.** Do a dedicated literature search for hidden-degree
   trapdoors combined with eigenbases. Compare with terSIDH and binSIDH
   (Basso–Fouotsa 2023), masked-degree SIDH variants, and the
   Castryck–Vercauteren backdoor discussion. None has been done beyond the
   surveys cited in the research note.
2. **Independent cryptanalysis** of §6.3, the partial knowledge of τ up to
   f, and of §6.4, the mask-guessing bound and Robert's N² > d condition, by
   someone other than the authors.
3. **A second implementation**, or at least an independent check of the
   known-answer vectors.
4. **Measurements owed:**
   * constant-time evidence at L3 and L5;
   * timings on a dedicated, pinned host;
   * a theta-coordinate chain for speed.
5. Re-run every number in §8 from its record before submission.
