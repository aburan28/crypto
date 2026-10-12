# ML-KEM: same-host reference, pq-crystals AVX2 (2026-10-05)

Stage: **reference measurement** (AGENTS.md §1). No code changes. This fixes
the boundary that `pqc::fast::ml_kem` (as merged in #1356) is measured
against, which `docs/pqc-speed.md` listed as "not done, not claimed".

## Reference

[pq-crystals/kyber](https://github.com/pq-crystals/kyber) at
`3edd5af5991927164edd4aacebfcbee00b8064e7` (the FIPS 203 version), `avx2/`
directory, its own `Makefile` and `test/test_speed{512,768,1024}`: hand-
written AVX2 assembly NTT, inverse NTT and pointwise product, XKCP four-way
Keccak, permuted ("shuffled") NTT coefficient order. Binary hashes in
`raw/pqcrystals_build.txt`.

Two properties of the build matter when reading the numbers:

- **It is compiled with `-march=native`.** On this host that lets gcc use
  AVX-512 in the C parts (notably XKCP's four-way Keccak, which gets
  `vprolq`). Ours chooses at run time from one portable binary. So the
  comparison is "their best build for this machine" against "our shipped
  binary", which is the right reference for a target but slightly favours
  them on the hashing.
- **Its assembly is AVX2 only.** The NTT, inverse NTT and pointwise product
  use no AVX-512 even here.

`kyber_keypair_derand` and `kyber_encaps_derand` take caller-supplied
randomness, which is what our `ml_kem_*_internal` functions are, so those
are the rows compared. Decapsulation is deterministic in both.

## Host

Same as `research/pqc_ml_kem_system_level_20261005/`: Intel Xeon @ 2.10 GHz
(family 6 model 207, KVM), 2 vCPUs, AVX2 + AVX-512 F/VL/BW/VBMI2, rustc
1.97.0, gcc 13.3.0. Every run through `tools/isolated_bench.py run --cpus 1`;
12 runs in `raw/target.jsonl`, 0 contended. Two earlier runs
(`raw/pqcrystals_first_runs.jsonl`) were contended and are not used.

## Method

Three rounds; each round runs pq-crystals at 512, 768 and 1024, then our
harness (`mbench cycles 3`, the same binary as the merged round's evidence).
Each cell below is the median over the three rounds of each tool's own
per-run median (pq-crystals: 1000 calls; ours: 600 samples). Per-round
values are in `raw/target_pqc.txt` and `raw/target_ours.txt`.

## The target

Kilocycles. "ratio" is ours ÷ pq-crystals: above 1 means we are slower.

| set | operation | pq-crystals AVX2 | ours | ratio | ours, prepared key | ratio |
|---|---|---:|---:|---:|---:|---:|
| 512 | keygen | 10.55 | 17.2 | 1.63 | 18.7 | 1.77 |
| 512 | encaps | 11.74 | 15.9 | 1.35 | 8.2 | 0.70 |
| 512 | decaps | 13.64 | 21.2 | 1.55 | 16.7 | 1.22 |
| 768 | keygen | 18.10 | 28.1 | 1.55 | 29.5 | 1.63 |
| 768 | encaps | 18.48 | 26.0 | 1.41 | 9.7 | 0.52 |
| 768 | decaps | 21.20 | 33.2 | 1.57 | 21.8 | 1.03 |
| 1024 | keygen | 23.93 | 36.1 | 1.51 | 37.5 | 1.57 |
| 1024 | encaps | 25.17 | 35.8 | 1.42 | 13.2 | 0.52 |
| 1024 | decaps | 29.19 | 44.4 | 1.52 | 29.2 | 1.00 |

Correctness of the reference was not re-checked here beyond its own build;
it is the reference implementation of FIPS 203. Ours is held byte-for-byte
to the ACVP-validated `pqc::ml_kem` by the merged tests.

Reading: on the operations both implement, **we are 1.35–1.63× behind**.
We are ahead only with a prepared key, which pq-crystals does not offer;
that is a different API, not a faster kernel, and pq-crystals would gain
from the same cache. Prepared decapsulation merely ties.

## Where the gap is: kernels

Cycles, median of 2000 calls (ours: `raw/kernel_cycles.rs`, run in the
scratch crate whose module is identical to the merged one) or 1000 calls
(pq-crystals' own `test_speed` rows), ML-KEM-768 unless noted.

| kernel | pq-crystals | ours | ratio |
|---|---:|---:|---:|
| NTT | 164 | 292 | 1.78 |
| inverse NTT | 174 | 328 | 1.89 |
| pointwise product, k = 3, with reduction | 184 | 358 | 1.95 |
| `gen_a` (expand `Â`), k = 2 | 2,530 | 3,188 | 1.26 |
| `gen_a`, k = 3 | 7,202 | 9,446 | 1.31 |
| `gen_a`, k = 4 | 10,078 | 12,620 | 1.25 |

Ours on the same run, for the rest of the budget: one η = 2 noise
polynomial 976; seven batched 3,238; compress d = 10 412; decompress
d = 10 326; `SHA3-256(ek)` at k = 3 5,634.

The ring arithmetic is ~1.8–1.95× slower per call and costs us ~3.6
kilocycles per ML-KEM-768 encapsulation (3 NTTs, 4 inverse NTTs, 4
pointwise products), ~1.7 more than the reference; `gen_a` is ~2.2 kilocycles behind. Together
they explain roughly half of the ~7.5 kilocycle encapsulation gap at
ML-KEM-768. The other half is not yet attributed; candidates are the
per-call allocation and zero-fill of buffers, lane extraction in the
four-way sponge, scalar CBD and packing, and the single sponges. It needs
a profile before it gets a fix.

## Success condition for the next round

Stated before any further change: the next round succeeds if, on this host
and against this table, ML-KEM-768 unprepared keygen, encaps and decaps
each reach a ratio ≤ 1.10 to pq-crystals AVX2, with byte-identical output.
It is abandoned for a given kernel if three attempts leave its ratio above
1.5. Prepared-key rows are reported but are not the success condition,
because the reference has no such API.

## Plan, in order of measured gap × share of time

1. **Attribute the unexplained half** of the gap (callgrind on the AVX2
   path for both, since valgrind cannot run the `-march=native` build:
   rebuild pq-crystals with `-mavx2 -mbmi2` for that profile only).
2. **`gen_a`:** lane extraction and the parse buffer, against XKCP's
   squeeze.
3. **Ring arithmetic:** stop restoring standard order between the short
   NTT layers (pq-crystals' permuted order, with pack/unpack at the byte
   boundary), lazy reduction, and keep the byte formats; this gives up
   "bit-identical to scalar" for "identical bytes out", so the tests move
   from array equality to encoded-output equality.
4. Then the prepared-key extras from the earlier list (cached Montgomery
   companions of `Â`), which pq-crystals cannot match without a cache.
