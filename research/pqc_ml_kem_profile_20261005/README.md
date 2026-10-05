# ML-KEM: where the gap to pq-crystals is (profile, 2026-10-05)

Stage: **stage diagnostic** (AGENTS.md §2). No code changes. This is step 1
of the plan in [`../pqc_ml_kem_reference_20261005/`](../pqc_ml_kem_reference_20261005/README.md):
attribute the gap between `pqc::fast::ml_kem` (main, after #1356) and
pq-crystals' AVX2 Kyber before changing anything.

## Method

Two independent views, because each misleads alone.

- **Instructions** (callgrind; deterministic, independent of the VM's
  neighbours). Valgrind does not emulate AVX-512, so both sides run their
  AVX2 path: ours by runtime dispatch, pq-crystals rebuilt with
  `-mavx2 -mbmi2 -mpopcnt` instead of its Makefile's `-march=native`. Per
  operation = (600-call run − 100-call run) / 500, so setup cancels.
  Drivers: `harness/pqc_driver.c` and the repository harness `mbench count`.
  Exclusive per-function tables: `raw/callgrind/*.no.agg`.
- **Cycles per phase** (`rdtsc`, isolated). Ours from
  `harness/kernel_cycles.rs` (per-call median of 2000); pq-crystals from
  its own `test_speed` rows, AVX2-only build. Both are same-key,
  back-to-back, per-call medians.

Every timed run went through `tools/isolated_bench.py run --cpus 1`; none
was marked contended.

## Result 1: the instruction gap is ~1.4–1.5×, and it is spread out

ML-KEM-768, instructions per operation, AVX2 path both sides:

| operation | ours | pq-crystals AVX2 | ratio |
|---|---:|---:|---:|
| keygen | 216,537 | 144,773 | 1.50 |
| encaps | 204,992 | 146,700 | 1.40 |
| decaps | 236,152 | 162,746 | 1.45 |

Encapsulation, by bucket (exclusive instructions; both sides make the same
number of four-way permutation calls, eight, and the same three scalar
`SHAKE128` blocks for the ninth matrix entry):

| bucket | ours | pq-crystals | Δ |
|---|---:|---:|---:|
| zero-fill and copies (`memset`, `memcpy`) | 23,569 | <800 | **+22.8k** |
| four-way Keccak permutation, 8 calls | 49,457 | 40,496 | **+9.0k** |
| pack, compress, decode `t̂` | 10,309 | ≈1,500 | **+8.8k** |
| noise: CBD and PRF glue | 6,553 | ≈1,000 | **+5.6k** |
| ring arithmetic (NTT, inverse, product, reduce) | 12,752 | 7,538 | **+5.2k** |
| rejection sampling and matrix glue | 12,312 | 7,509 | **+4.8k** |
| scalar Keccak (`H(ek)`, `G`, 9th entry) | 80,795 | 77,194 | +3.6k |
| total | 204,992 | 146,700 | +58.3k |

Two cautions on reading this table as time:

- valgrind counts each byte of `rep stosb` as an instruction, so the
  zero-fill row overstates its cycle cost; zeroing the `KMAX²`-polynomial
  matrix costs 120 cycles measured, not thousands.
- our four-way Keccak is ~22% more instructions per call than XKCP's on the
  AVX2 path; on this AVX-512 host the shipped code uses the AVX-512 kernel
  instead, which is why our matrix expansion is *faster* than
  pq-crystals' AVX2-only build in cycles (below).

Decapsulation adds one more large item: the ciphertext comparison and
implicit-rejection select through `subtle` (`ml_kem_decaps` and `<u8>`
exclusive, ~13k instructions), which compares byte by byte; pq-crystals'
`verify.c` does it 32 bytes at a time.

## Result 2: in cycles, the gap is the small kernels

ML-KEM-768 phases, cycles per call:

| phase | calls per encaps | ours | pq-crystals AVX2 | Δ per encaps |
|---|---:|---:|---:|---:|
| compress + pack `u` (d = 10) | 3 | 358 | 69 (208 for three) | **+866** |
| inverse NTT | 4 | 280 | 156 | **+496** |
| pointwise product, k = 3 | 4 | 320 | 164 | **+624** |
| compress + pack `v` (d = 4) | 1 | 338 | 32 | **+306** |
| NTT | 3 | 240 | 144 | +288 |
| `Decompress_1(m)` | 1 | 270 | 30 | +240 |
| matrix expansion, k = 3 | 1 | 8,290 | 9,290 | **−1,000** |

Decapsulation also runs: decompress + unpack `u` 3 × 290 against 78 for
three (**+792**), decompress `v` 154 against 28, `Compress_1` (to message)
330 against 30 (**+300**), decode `ŝ` 690, and everything in the encaps
column again for re-encryption.

The pattern: our ring arithmetic is ~1.7–2.0× theirs per call, and our
byte-boundary kernels (pack, unpack, compress, decompress, message
encode/decode) are **5–11×** theirs. pq-crystals does those in AVX2 end to
end; ours computes in vectors and then packs in scalar code, or is scalar
throughout. Those kernels alone are ~1.4k cycles of an encapsulation and ~2.6k
of a decapsulation, before counting the decoding of `t̂` and `ŝ`.

## Result 3: the reference ratio has wide error bars on this host

Five interleaved rounds (pq-crystals native, pq-crystals AVX2-only, ours;
`raw/round2.tsv`), all uncontended by the tool's test, gave ML-KEM-768
encapsulation ratios, ours ÷ reference, per round:

| round | vs pq-crystals native | vs pq-crystals AVX2-only |
|---|---:|---:|
| 1 | 1.63 | 1.26 |
| 2 | 1.54 | 1.19 |
| 3 | 1.37 | 1.03 |
| 4 | 1.39 | 1.07 |
| 5 | 2.25 | 0.95 |

And the same code measured by two builds of the same harness differed by
~12% (`mbench` built for #1356: 26.0k; rebuilt today: 23.0k), the
build-to-build layout effect `docs/pqc-speed.md` already warns about at
small scale. A throughput loop (`harness/tput.rs`, `harness/pqc_tput.c`,
`raw/tput.tsv`) was noisier still.

So the 1.41× in #1376's table is one sample from a 1.37–1.63 spread against
the native build. The stable statement is the instruction ratio above
(1.40× on encapsulation, same AVX2 path on both sides), and against an
AVX2-only reference the cycle gap is smaller, roughly 1.0–1.3×. Class:
**accounting**. It does not change the plan, but the success condition
should be judged on instructions plus interleaved same-round cycle ratios,
not on one session's table.

## What to do, ranked by measured gap

1. **AVX2 pack/compress and unpack/decompress, fused**, including the
   one-bit message codecs and `t̂`/`ŝ` decoding. Largest cycle gap per
   call (5–11×); ~1.4k cycles per encapsulation, ~2.6k per decapsulation,
   plus key decoding.
2. **Stop zero-filling and copying polynomials**: sample straight into the
   matrix and noise vectors, no zero-initialised scratch, no `[Poly; 4]`
   returns. ~23k instructions; cycle value to be measured.
3. **Vectorised constant-time compare and select** in decapsulation,
   still branch-free. ~13k instructions.
4. **Ring arithmetic**: fewer shuffles between the short NTT layers, lazy
   reduction. ~1.4k cycles per encapsulation.
5. **Vector CBD** (~4.6k instructions) and **register pressure in the AVX2
   Keccak** (~22% more instructions than XKCP; AVX-512 hosts unaffected).

Items 1–3 are bookkeeping, not arithmetic, and together are most of the
gap. That inverts the expectation the earlier rounds started from.
