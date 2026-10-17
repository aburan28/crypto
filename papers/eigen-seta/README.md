# Eigen-Séta: an isogeny KEM whose trapdoor is a hidden-degree endomorphism

`eigen_seta.tex` specifies Eigen-Séta v2, proves that it is correct, derives
parameters for λ = 128, 192 and 256, analyses the attacks we know, and
describes and measures the Rust implementation in `crypto:eigen_seta/`.

**Status.** Research draft, revision 2 (2026-10-17; first draft 2026-10-11),
anonymized, not peer reviewed, not submitted and not posted. Security is a
conjecture supported only by the analysis in the paper. Restore the authors
with `make deanonymize` after creating an `AUTHORS` file (see
`../AUTHORS.example`).

**What changed in revision 2.** An assessment of the scheme against the
maintained implementation, the sibling repository `eigenx`, was added as §9,
along with:

* the public-key degree test (Lemma 8, §6.3): one candidate f is tested with a
  Kani chain of about 2·log₂ f steps on E_A × E_A, with no ciphertext. It is
  measured in §8.4: the true degree was accepted and all 158 wrong
  candidates were rejected at L1, L3 and L5. Its arithmetic core is checked
  in Lean;
* the NIST-category conjecture of revision 1 is **withdrawn in its unqualified
  form** and replaced by two tiers (§6.5). Tier A is gate count under
  MAXDEPTH: L1, L3 and L5 conjectured to meet categories 1, 3 and 5. Tier B
  is at least 2^λ quantum queries: only L5 (as level 1), Q3 and Q5 meet it;
* the Q3 and Q5 parameter sets, NUMS primes, a cost table for all five sets;
* the hidden orientation of E_A named as the largest structural risk;
* the newer literature: Wesolowski p^{1/3}, Udovenko p^{2/5}, the
  one-endomorphism result of Herlédan Le Merdy and Wesolowski, self-pairings,
  adaptive attacks on FESTA, faults;
* pointers to the maintained code, which decapsulates with a theta chain.
  The revision-1 measurements in §8.1–8.3 are kept and labelled historical
  (Mumford build).

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
python figures/plot_assessment.py ../../research/eigen_seta_assessment_20261017 \
    ../../research/eigen_seta_assessment_20261017/figures   # then copy pk_test.pdf here
```

## Reproduction

| What | Where |
|---|---|
| maintained code, tests, KATs, cost model | `eigenx:` at commit `94fe2d9` (`cargo test --release`, `cargo run --release --example security_estimate`) |
| public-key degree test (§6.3, §8.4) | `eigenx:examples/pk_degree_test.rs`; records `crypto:research/eigen_seta_assessment_20261017/records/run-20261017T1623Z-pk-degree-test/` |
| Lean check of Lemma 8's arithmetic | `crypto:research/eigen_seta_assessment_20261017/lean/` (`tools/lean/proof-gate.sh <file>`) |
| assessment report (graphs, PDF) | `crypto:research/eigen_seta_assessment_20261017/REPORT.md` |
| historical code (revision 1) | `crypto:eigen_seta/` |
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
   * constant-time evidence on the theta build, and at L3 and L5;
   * timings on a dedicated, pinned host (the revision-2 timings are on a
     shared 4-core VM);
   * a Q5 run on the reference host.
5. **Hidden orientation.** Instantiate the level-structure frameworks
   (De Feo–Fouotsa–Panny 2024/459, Arpin et al. 2024/1172) and sesquilinear
   pairings against the public key (§9, risk 1).
6. Re-run every number in §8 from its record before submission.
