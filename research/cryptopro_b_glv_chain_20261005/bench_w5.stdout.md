# cryptopro-b-glv-bench

curve: GOST R 34.10-2001 CryptoPro-B (RFC 4357 section 11.4.4), p = 2^255 + 3225, prime order n, cofactor 1
arithmetic: variable-time Jacobian a = -3 (dbl-2001-b, add-2007-bl, madd-2007-bl), Montgomery-form field, shared by every arm
protocol: M = 256 seeded (scalar, point) pairs, R = 20 rounds, w = 5, seed = 20261005, arm order alternates between rounds; timed output is a Jacobian point (no final inversion in any arm)
host: Intel(R) Xeon(R) Processor @ 2.10GHz; flags: pclmulqdq sse4_2 popcnt aes avx bmi1 avx2 bmi2 avx512f adx; logical cpus: 1; kernel: 6.18.44-fc-v70; linux x86_64
build: rustc 1.97.0 (2d8144b78 2026-07-07); profile release; commit e1df564168943f7d03f9462503feb63b218fa392 (dirty)
timing wall time: 2.2 s

| arm | class | ns per scalar multiplication, median of 20 rounds | min of 20 rounds | ratio baseline/arm (median) | ratio baseline/arm (min) | correctness |
|:--|:--|--:|--:|--:|--:|:--|
| baseline: width-5 NAF (affine odd-multiple table, mixed additions) | reference | 111158 | 107220 | 1.000 | 1.000 | baseline == textbook BigUint reference on all 256 pairs: yes |
| GLV-2 total via chain 5·5·7 (element 4 + 1 w, degree 175): decomposition + phi(P) + interleaved width-5 NAF | engineering | 88759 | 80199 | 1.252 | 1.337 | GLV == baseline on all 256 pairs: yes |
| stage diagnostic: decomposition alone, basis of chain 5·5·7 (num-bigint Babai rounding) | stage diagnostic | 1006 | 875 | — (stage) | — (stage) | k1 + k2·λ ≡ k (mod n) and |k1|,|k2| ≤ 128 bits on all 256 scalars: yes |
| stage diagnostic: phi(P) alone, chain 5·5·7 (element 4 + 1 w, degree 175; projective Horner, no inversion) | stage diagnostic | 4269 | 3841 | — (stage) | — (stage) | phi(P) == λ·P on all 256 points: yes |
| GLV-2 total via chain 5·31 (element 0 + 1 w, degree 155): decomposition + phi(P) + interleaved width-5 NAF | engineering | 91799 | 84544 | 1.211 | 1.268 | GLV == baseline on all 256 pairs: yes |
| stage diagnostic: decomposition alone, basis of chain 5·31 (num-bigint Babai rounding) | stage diagnostic | 1007 | 863 | — (stage) | — (stage) | k1 + k2·λ ≡ k (mod n) and |k1|,|k2| ≤ 128 bits on all 256 scalars: yes |
| stage diagnostic: phi(P) alone, chain 5·31 (element 0 + 1 w, degree 155; projective Horner, no inversion) | stage diagnostic | 9510 | 8479 | — (stage) | — (stage) | phi(P) == λ·P on all 256 points: yes |
| A/A control: the baseline arm timed again as a separate arm (identical code and inputs; its ratio to the baseline is the noise floor) | A/A control | 118769 | 108890 | 0.936 | 0.985 | same routine as the baseline, so baseline == textbook BigUint reference on all 256 pairs: yes |

all correctness checks passed: yes
json record written to research/cryptopro_b_glv_chain_20261005/bench_w5.json
