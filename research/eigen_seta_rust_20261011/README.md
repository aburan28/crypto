# Eigen-Séta v2 in Rust: KEM and authenticated PKE at L1, L3 and L5 sizes

Results report, 2026-10-11. The PDF is `report.pdf`. The code is in the
crate `eigen_seta/`, whose README describes the scheme and whose
`SIDE_CHANNELS.md` holds the full audit. This report gathers the evidence
for what the accompanying PR claims and states the limits of each claim.

## Claims and their status

| # | Claim | Status | Evidence |
|---|---|---|---|
| 1 | The Rust KEM decapsulates correctly at Toy16, Toy32, L1, L3 and L5. | measured: 15/15 trials in the benchmark run, plus earlier development runs | `reports/run-20261011T1813Z-bench/` |
| 2 | Public key and ciphertext are 1344, 1992 and 2632 bytes at L1, L3 and L5. Decapsulation takes 13–15 s, 52–54 s and 120–124 s on one core of the recorded host. | measured, one host, n = 3, 2 and 2 | same |
| 3 | The original secret phase of decapsulation ran a key-independent and ciphertext-independent instruction stream, apart from glibc's allocator. | measured at Toy16 for 4 keys × 3 ciphertexts with callgrind on commit `19c641d12`; **not rerun after the later scalar-recovery and input-validation edits**. This is historical evidence, not a proof for the current build. | `reports/run-20261011T1827Z-sidechannel/` |
| 4 | No timing difference is detected between a fixed ciphertext and fresh ones, including FO-rejected ones, for one key. | measured with a Welch t-test, 1,000 samples per class: \|t\| ≤ 0.84, against a threshold of 4.5 | same |
| 5 | The parameter sets target λ = 128, 192 and 256 classical bits. | **conjectured**, with heuristic arguments in `eigen_seta/README.md` §3. No independent cryptanalysis exists. Quantum security is about λ/2, as for FESTA. | — |
| 6 | The toy `walk` variant in the merged prototype is insecure. | measured: an outsider decrypted 5/5 ciphertexts from public data | correction notices in `RESEARCH_ENDRING_PKE_STRATEGY.md` and the prototype README |
| 7 | The KEM–DEM PKE API encrypts and decrypts a 39-byte plaintext with associated data at L1, L3 and L5. | measured: 1/1 round trip per size on the later build; Toy16 tamper/context/truncation tests pass. This is functionality, not security validation. | `reports/run-20261011T1956Z-pke-l1/`, `reports/run-20261011T1958Z-pke-l3-l5/`, `eigen_seta/tests/pke.rs` |

Two defects found during this work are fixed in the same PR:

* **An LLVM-introduced conditional in safegcd.** The masked add in
  `normalize` was compiled into a conditional skip, so the timing depended
  on the sign of a secret-dependent intermediate. The callgrind probe found
  it, and `black_box` on the masks fixed it. The counts above were taken
  after the fix.
* **A key-generation hang.** Tonelli–Shanks accepted z = 0 as a non-residue.
  A regression test covers it.
* **Degenerate scalar recovery.** The Okeya–Sakurai formula has a zero
  denominator when `[k+1]P = O`. The scalar multiplication now selects
  `−P` in that case, with tests on 2-power and 3-power torsion. The compiled
  side-channel probe has not been repeated on this later build.
* **Unbounded basis search on malformed public curves.** The public search
  is now bounded, and encapsulation validates the public 2-power basis.

## Results

### Correctness, sizes and speed

The benchmark run is `run-20261011T1813Z-bench`, on an Intel Xeon at
2.80 GHz in a shared VM, built with rustc 1.97.0. Every trial's decapsulated
key matched the encapsulated one.

| Set | p bits | pk = ct (B) | keygen (s), n = 1 | encaps (s) | decaps (s) | correct |
|---|---|---|---|---|---|---|
| Toy16 | 216 | 216 | 0.38 | 0.005–0.008 | 0.12–0.19 | 5/5 |
| Toy32 | 369 | 376 | 2.3 | 0.009–0.012 | 0.45–0.59 | 3/3 |
| L1 | 1341 | 1344 | 32.3 | 0.24–0.26 | 13.4–15.0 | 3/3 |
| L3 | 1991 | 1992 | 118.6 | 1.03 | 51.8–53.9 | 2/2 |
| L5 | 2630 | 2632 | 271.8 | 2.17–2.20 | 120.2–124.0 | 2/2 |

![KEM wall time](figures/kem_times.svg)

Measured in the same binary and run, safegcd inversion is 3.9× to 140×
faster than Fermat inversion, from Toy16 to L5. Field multiplication scales
as expected. No pre-change baseline binary was kept for the KEM-level
timings, so this report makes **no KEM-level speedup claim**. AGENTS.md §10
requires a recorded baseline for one.

![Field operations](figures/field_ops.svg)

### Follow-up: authenticated plaintext encryption

The new `pke::EigenSetaPke` combines a fresh Eigen-Séta KEM encapsulation
with ChaCha20-Poly1305. The complete ciphertext contains the KEM ciphertext,
a 12-byte nonce, the encrypted message and a 16-byte tag. An integration test
checks empty and nonempty plaintext, changed associated data, tampered KEM
component/nonce/tag, truncation, and invalid public keys at Toy16. The
command-line example completed one 39-byte round trip at each full-size set.

| Set | KeyGen (s) | Encrypt (s) | Decrypt (s) | pk (B) | complete PKE ct (B) | correct |
|---|---:|---:|---:|---:|---:|---:|
| L1 | 25.044 | 0.3897 | 11.208 | 1344 | 1411 | 1/1 |
| L3 | 186.427 | 1.1051 | 47.138 | 1992 | 2059 | 1/1 |
| L5 | 298.906 | 2.2990 | 104.223 | 2632 | 2699 | 1/1 |

Each time is one measured wall-clock sample on an AMD EPYC 7763 shared VM
with rustc 1.99.0; no spread was recorded. L3/L5 were concurrent and other
validation was active. These are **descriptive correctness-run timings**, not
an isolated benchmark or a speed comparison with the older Intel Xeon run.
Source records: `reports/run-20261011T1956Z-pke-l1/records/run_output.txt`
and `reports/run-20261011T1958Z-pke-l3-l5/records/run_output.txt`.

![PKE round-trip wall time](figures/pke_roundtrips.svg)

The current build also completed new KEM round trips at L1, L3 and L5
(`reports/run-20261011T1944Z-l1-fix/`,
`reports/run-20261011T1946Z-l3-l5-fix/`, and
`reports/run-20261011T1955Z-final-l1/`). The scalar-recovery and
public-input edits are covered by targeted tests. The original random-curve
two-halves construction remains distinct and unimplemented; these results
belong to the secret-degree single-chain v2 variant. Security is conjectural.

### Constant-time evidence

`run-20261011T1827Z-sidechannel` covers the callgrind instruction counts of
`decaps_core` on the earlier commit `19c641d12`, at Toy16. Each of 4 key seeds has 2 valid
ciphertexts and 1 foreign, FO-rejected one. The table of values is in
`SIDE_CHANNELS.md` §5 and in the run report. The non-allocator count was
identical in all 12 runs. The glibc allocator count varies with the heap
state left by the earlier variable-time phases.

![Constant-time probe](figures/ct_probe.svg)

Limits:

* Instruction counts do not cover memory-address leakage. The code has no
  secret-indexed tables, which is an argument, not a measurement.
* The check was run at Toy16 only. The code is generic over the limb count
  N, so the same control flow is expected at L1–L5, but that was not
  measured.
* The result covers this compiler and target only.
* It does not cover the later scalar-recovery fix or the PKE layer; no
  compiled-code constant-time claim is made for this PR's final build.
* Key generation is variable time by design. See `SIDE_CHANNELS.md`
  residual issue R1.

## Graphs checked and left unchanged

These canonical graphs were checked:

* `docs/index-calculus-scoreboard.html`
* `docs/ic/progress-timeline.json`
* `docs/performance-gains/`
* `docs/curves/`
* `figures/`

None of them covers isogeny-based KEMs or this crate. `docs/curves/`
mentions isogenies only as curve-transfer links in the discrete-log curve
registry, and this work adds no registry curve. So none of them changed. This
report's graphs live in `figures/`. The historical plots were made from run
records by `figures/plot.py`; the follow-up PKE SVG is drawn from
`figures/pke_roundtrips.csv`, whose values are copied from its two run records.

## Files

* `README.md`, `report.typ`, `report.pdf`: this report.
* `figures/`: `kem_times.svg`, `field_ops.svg`, `ct_probe.svg`,
  `pke_roundtrips.svg` with its CSV data, and the earlier plotting script.
* `reports/run-20261011T1813Z-bench/`: the benchmark run report, its
  records and its PDF.
* `reports/run-20261011T1827Z-sidechannel/`: the side-channel run report,
  its records (callgrind outputs xz-compressed, probe outputs, dudect
  output) and its PDF.
* `reports/run-20261011T1944Z-l1-fix/`, `run-20261011T1946Z-l3-l5-fix/`,
  `run-20261011T1955Z-final-l1/`: follow-up KEM run reports, records,
  charts and PDFs.
* `reports/run-20261011T1956Z-pke-l1/` and
  `run-20261011T1958Z-pke-l3-l5/`: PKE run reports, records, charts and PDFs.
