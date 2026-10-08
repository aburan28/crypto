# ML-KEM: system-level optimisation round (2026-10-05)

Stage: **engineering** (AGENTS.md §3). This changes the cost of a correct
implementation of a fixed algorithm; it is not a cryptanalytic result and no
boundary moved. Every figure here is kilocycles (`rdtsc`) or instructions per
operation, and every row is byte-identical to the ACVP-validated reference.

## Question and protocol

Hypothesis, from Sayed, Taha and Nijjer, *System-Level Optimization Beyond
Cryptographic Kernels: An ML-KEM Case Study on Arm Cortex-M7*
(arXiv:2610.01960): once arithmetic kernels are tuned, the largest remaining
cost in ML-KEM is (a) public-only work repeated per call, chiefly expanding
`Â` and hashing `ek`, and (b) the surrounding execution system. On a
Cortex-M7 they measure up to 74.6% / 58.8% fewer encaps / decaps cycles from
public-data reuse and at most ~2.5% from everything else.

The question for this repository: does the same hold for
`pqc::fast::ml_kem` on x86-64, and what do the x86 counterparts of the
paper's "execution system" knobs (wide vector units in place of TCM and
peripherals) buy on top?

- **Reference (unmodified baseline).** `main` at `9f5afa23`, built and kept
  as a binary before any change (`raw/binaries.sha256`).
- **Frozen inputs.** 16 deterministic `(d, z, m)` seed triples per parameter
  set from the xorshift generator in `harness/mbench.rs` (tags 1, 2, 3),
  rotated across samples so the rejection samplers do not run one path.
- **Correctness gate.** Before any timing, every side must produce identical
  keys, ciphertexts and shared secrets on every frozen input (asserted in the
  harness fixture). In the repository harness, the prepared ciphertext must
  decapsulate under the *reference* implementation to the same secret.
- **Stop condition.** Any byte difference against the reference, or any
  hardware class (AVX-512, AVX2-only, scalar) slower than baseline.
- **Cost accounting.** Prepared-key rows exclude the one-time preparation,
  which is reported as its own row (`prep-ek`, `prep-dk`). Key generation
  that captures the cache is reported in full.

Honesty note: this protocol was written up during the work rather than
committed before the first run. The baseline binary and its A/A spread were
recorded before the first candidate measurement, and the inputs above were
fixed then.

## Host

| | |
|---|---|
| CPU | Intel Xeon @ 2.10 GHz (family 6 model 207, KVM guest), 2 vCPUs |
| features | AVX2, AVX-512 F/VL/BW/VBMI2, BMI2, AES-NI, SHA-NI |
| memory | 7 GiB |
| OS | Linux 6.18 x86_64 |
| toolchain | rustc 1.97.0, `--release`, no `target-cpu` flags |
| isolation | `tools/isolated_bench.py run --cpus 1`, every timed run; 0 contended |

## Commands

```sh
# scratch crate depending on this repository by path; side a = main's
# pqc::fast::ml_kem, b = candidate, c = candidate with prepared keys
python3 tools/isolated_bench.py run --cpus 1 --out raw/stepN.jsonl -- \
    mbench cycles 7                 # ABC interleaved, 200 samples per block
valgrind --tool=callgrind mbench count <a|b|c> <set> <op> <100|600>
                                    # difference of the two = 500 operations
python3 tools/isolated_bench.py run --cpus 1 --out raw/pqc_speed_ab.jsonl -- \
    pqc_speed ml-kem                # base, new, base, new, ... three rounds
```

`harness/` holds the scratch crate's manifest and harness source. The
candidate side was compiled from the scratch copy of the module whose final
state is this PR's `src/pqc/fast`; the only differences are module paths and
two lint fixes.

## A/A

Baseline against an identical copy, interleaved: medians within ±2% on every
row (`raw/aa.jsonl` conditions; table in `tables.md`). Differences below 3%
are not claimed.

## Result

See `tables.md` for every run. Final interleaved table, ML-KEM-768,
kilocycles, median of 1400 samples per cell:

| operation | main | this PR | prepared key | main/PR | main/prepared |
|---|---:|---:|---:|---:|---:|
| keygen | 83.1 | 28.0 | 29.6 | 2.97× | 2.81× |
| encaps | 96.1 | 25.9 | 9.8 | 3.72× | 9.81× |
| decaps | 119.8 | 33.1 | 22.3 | 3.62× | 5.38× |

Instructions per operation (callgrind; valgrind does not emulate AVX-512, so
this is the AVX2 code path): ML-KEM-768 encaps 527,337 → 204,964 →
56,993 prepared.

Decision: merge. The paper's ordering holds on x86 for the *first* step
(public-data reuse is the largest single lever, 2.18× on encaps), but unlike
the M7, the execution-system levers are not marginal here: SIMD Keccak and
SIMD ring arithmetic are worth more than the cache once the cache exists.
