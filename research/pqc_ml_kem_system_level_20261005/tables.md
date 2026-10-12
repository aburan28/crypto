# Raw tables, in the order they were run

Kilocycles (`rdtsc`, reference clock), median per cell. `base` is `main`
(`9f5afa23`) `pqc::fast::ml_kem`; `cand` is the candidate's unprepared API;
`prepared` is the candidate with a prepared key. Sides run interleaved
(a block of 200 samples of each, repeated per round) in one process, pinned
by `tools/isolated_bench.py`; the matching condition records are in `raw/`.
Superseded rows are kept: each step's table is the state after that step.

## A/A (candidate identical to base)

Not isolated (first run, before the tool was set up); medians agree within
±2%, which is taken as the noise floor.

| set | op | base med | cand med | ratio |
|---|---|---:|---:|---:|
| 512 | keygen | 45.7 | 46.5 | 0.98 |
| 512 | encaps | 53.0 | 52.5 | 1.01 |
| 512 | decaps | 71.8 | 70.6 | 1.02 |
| 768 | keygen | 82.0 | 83.2 | 0.98 |
| 768 | encaps | 92.6 | 93.5 | 0.99 |
| 768 | decaps | 117.6 | 118.9 | 0.99 |
| 1024 | keygen | 133.0 | 133.4 | 1.00 |
| 1024 | encaps | 146.9 | 147.8 | 0.99 |
| 1024 | decaps | 178.3 | 177.8 | 1.00 |

A second, isolated A/A (the first record in `raw/step1.jsonl`, label
`prepared`: a stale harness binary that still compared two identical copies
of the baseline, run by mistake) gave 0.97–1.01 on every row.

## Step 1: public-data reuse (`raw/step1.jsonl`)

| set | op | base | cand | prepared | base/cand | base/prepared |
|---|---|---:|---:|---:|---:|---:|
| 512 | keygen | 48.8 | 48.0 | 47.4 | 1.02 | 1.03 |
| 512 | encaps | 55.2 | 53.2 | 30.9 | 1.04 | 1.79 |
| 512 | decaps | 71.8 | 70.7 | 50.8 | 1.02 | 1.41 |
| 768 | keygen | 84.4 | 82.8 | 82.1 | 1.02 | 1.03 |
| 768 | encaps | 94.4 | 93.6 | 43.4 | 1.01 | 2.18 |
| 768 | decaps | 117.7 | 117.6 | 70.7 | 1.00 | 1.66 |
| 1024 | keygen | 134.5 | 132.8 | 132.8 | 1.01 | 1.01 |
| 1024 | encaps | 146.7 | 145.7 | 59.7 | 1.01 | 2.46 |
| 1024 | decaps | 178.3 | 177.0 | 95.5 | 1.01 | 1.87 |

Preparation: ek 20.4 / 47.7 / 84.0, dk 18.8 / 43.4 / 78.6 (512/768/1024).
Key generation with capture costs nothing measurable; on the M7 the paper
measured +2.5–5.2%, because there the capture is a serialisation.

## Step 2: four-way Keccak (`raw/step2.jsonl`)

| set | op | base | cand | prepared | base/cand | base/prepared |
|---|---|---:|---:|---:|---:|---:|
| 512 | keygen | 48.5 | 33.5 | 35.2 | 1.45 | 1.38 |
| 512 | encaps | 55.1 | 39.2 | 27.0 | 1.40 | 2.04 |
| 512 | decaps | 73.1 | 56.3 | 47.1 | 1.30 | 1.55 |
| 768 | keygen | 84.0 | 56.4 | 57.8 | 1.49 | 1.45 |
| 768 | encaps | 95.3 | 62.8 | 38.7 | 1.52 | 2.46 |
| 768 | decaps | 118.0 | 85.8 | 66.2 | 1.37 | 1.78 |
| 1024 | keygen | 134.5 | 77.2 | 78.6 | 1.74 | 1.71 |
| 1024 | encaps | 147.9 | 86.2 | 54.2 | 1.72 | 2.73 |
| 1024 | decaps | 178.3 | 116.2 | 89.9 | 1.53 | 1.98 |

Permutation microbenchmark, cycles: scalar one state 882–949; AVX2 four
states 1504–1693; AVX-512VL four states 609–616.

## Step 3: AVX2 NTT, inverse NTT, pointwise product (`raw/step3.jsonl`)

| set | op | base | cand | prepared | base/cand | base/prepared |
|---|---|---:|---:|---:|---:|---:|
| 512 | keygen | 47.5 | 23.0 | 24.8 | 2.06 | 1.91 |
| 512 | encaps | 54.4 | 23.9 | 11.8 | 2.28 | 4.59 |
| 512 | decaps | 72.4 | 34.7 | 24.9 | 2.09 | 2.91 |
| 768 | keygen | 83.4 | 38.4 | 40.1 | 2.17 | 2.08 |
| 768 | encaps | 94.4 | 38.4 | 14.3 | 2.46 | 6.60 |
| 768 | decaps | 117.7 | 52.9 | 32.3 | 2.23 | 3.64 |
| 1024 | keygen | 133.2 | 49.3 | 51.2 | 2.70 | 2.60 |
| 1024 | encaps | 146.5 | 50.8 | 18.7 | 2.88 | 7.82 |
| 1024 | decaps | 179.1 | 70.3 | 42.5 | 2.55 | 4.21 |

## Step 4: grouped packing, single sponges in one AVX-512 lane (`raw/step4.jsonl`)

| set | op | base | cand | prepared | base/cand | base/prepared |
|---|---|---:|---:|---:|---:|---:|
| 512 | keygen | 47.9 | 17.7 | 19.3 | 2.71 | 2.48 |
| 512 | encaps | 55.2 | 18.1 | 9.7 | 3.05 | 5.71 |
| 512 | decaps | 73.6 | 24.1 | 18.4 | 3.05 | 4.00 |
| 768 | keygen | 83.6 | 30.4 | 31.9 | 2.75 | 2.62 |
| 768 | encaps | 94.3 | 30.3 | 11.8 | 3.11 | 7.97 |
| 768 | decaps | 118.6 | 38.5 | 24.1 | 3.08 | 4.92 |
| 1024 | keygen | 136.1 | 39.0 | 40.4 | 3.49 | 3.37 |
| 1024 | encaps | 148.1 | 41.1 | 16.4 | 3.60 | 9.05 |
| 1024 | decaps | 180.7 | 52.4 | 33.0 | 3.45 | 5.48 |

Single permutation: scalar 882–933 cycles; one lane of the AVX-512 kernel
631–714.

## Step 5: vector compression (`raw/step5.jsonl`, second run; the first is kept there and had a 141.2 base outlier on 768 encaps)

| set | op | base | cand | prepared | base/cand | base/prepared |
|---|---|---:|---:|---:|---:|---:|
| 512 | keygen | 46.0 | 18.5 | 19.6 | 2.50 | 2.35 |
| 512 | encaps | 51.8 | 16.6 | 8.2 | 3.12 | 6.32 |
| 512 | decaps | 68.8 | 21.8 | 16.8 | 3.16 | 4.11 |
| 768 | keygen | 81.7 | 30.8 | 32.1 | 2.65 | 2.55 |
| 768 | encaps | 91.9 | 28.1 | 9.8 | 3.27 | 9.40 |
| 768 | decaps | 116.4 | 35.2 | 21.8 | 3.31 | 5.34 |
| 1024 | keygen | 133.4 | 39.3 | 40.7 | 3.40 | 3.28 |
| 1024 | encaps | 146.2 | 37.3 | 13.0 | 3.92 | 11.28 |
| 1024 | decaps | 177.1 | 47.2 | 29.3 | 3.76 | 6.04 |

## Step 6: vector rejection sampling (`raw/step6.jsonl`) — final

| set | op | base | cand | prepared | base/cand | base/prepared |
|---|---|---:|---:|---:|---:|---:|
| 512 | keygen | 46.7 | 17.1 | 18.7 | 2.73 | 2.50 |
| 512 | encaps | 54.2 | 15.9 | 8.3 | 3.40 | 6.55 |
| 512 | decaps | 72.0 | 21.3 | 16.9 | 3.38 | 4.25 |
| 768 | keygen | 83.1 | 28.0 | 29.6 | 2.97 | 2.81 |
| 768 | encaps | 96.1 | 25.9 | 9.8 | 3.72 | 9.81 |
| 768 | decaps | 119.8 | 33.1 | 22.3 | 3.62 | 5.38 |
| 1024 | keygen | 136.3 | 36.5 | 37.8 | 3.74 | 3.61 |
| 1024 | encaps | 149.0 | 34.5 | 13.0 | 4.32 | 11.43 |
| 1024 | decaps | 180.5 | 44.4 | 30.0 | 4.07 | 6.01 |

Preparation: ek 7.7 / 16.2 / 21.6, dk 5.6 / 12.5 / 16.1.

## Per hardware class (`raw/classes.jsonl`)

The final code built with the backend pinned, interleaved against the same
baseline. The pins were a scratch-only build feature; the shipped code
chooses at run time.

AVX2 only (no AVX-512):

| set | op | base | cand | prepared | base/cand | base/prepared |
|---|---|---:|---:|---:|---:|---:|
| 512 | keygen | 53.5 | 23.3 | 27.7 | 2.30 | 1.93 |
| 512 | encaps | 54.3 | 22.6 | 10.3 | 2.40 | 5.28 |
| 512 | decaps | 71.1 | 28.0 | 20.4 | 2.54 | 3.49 |
| 768 | keygen | 82.9 | 36.9 | 40.5 | 2.25 | 2.04 |
| 768 | encaps | 91.6 | 36.5 | 11.7 | 2.51 | 7.84 |
| 768 | decaps | 116.3 | 43.7 | 26.2 | 2.66 | 4.44 |
| 1024 | keygen | 136.2 | 52.5 | 52.5 | 2.59 | 2.59 |
| 1024 | encaps | 145.7 | 52.0 | 14.9 | 2.80 | 9.80 |
| 1024 | decaps | 176.6 | 61.7 | 34.6 | 2.86 | 5.11 |

Scalar only (the path Arm64 and pre-AVX2 x86 take; measured here on x86):

| set | op | base | cand | prepared | base/cand | base/prepared |
|---|---|---:|---:|---:|---:|---:|
| 512 | keygen | 45.5 | 43.8 | 45.1 | 1.04 | 1.01 |
| 512 | encaps | 52.4 | 48.3 | 29.0 | 1.08 | 1.80 |
| 512 | decaps | 69.3 | 60.7 | 46.8 | 1.14 | 1.48 |
| 768 | keygen | 82.7 | 71.1 | 72.3 | 1.16 | 1.14 |
| 768 | encaps | 92.9 | 79.6 | 41.1 | 1.17 | 2.26 |
| 768 | decaps | 116.7 | 96.8 | 66.0 | 1.21 | 1.77 |
| 1024 | keygen | 138.4 | 114.4 | 117.8 | 1.21 | 1.17 |
| 1024 | encaps | 149.6 | 121.7 | 57.4 | 1.23 | 2.61 |
| 1024 | decaps | 178.1 | 144.7 | 90.2 | 1.23 | 1.97 |

## Instructions per operation (`raw/callgrind_ab.tsv`)

valgrind 3.x does not expose AVX-512, so these are the AVX2 path's counts.
Difference of a 600- and a 100-operation run, divided by 500, so the fixture
cancels.

| set | op | base | cand | prepared | base/cand | base/prepared |
|---|---|---:|---:|---:|---:|---:|
| 512 | keygen | 296,094 | 136,026 | 142,809 | 2.18 | 2.07 |
| 512 | encaps | 323,580 | 132,343 | 54,738 | 2.45 | 5.91 |
| 512 | decaps | 409,835 | 155,297 | 113,570 | 2.64 | 3.61 |
| 768 | keygen | 480,454 | 216,577 | 223,596 | 2.22 | 2.15 |
| 768 | encaps | 527,337 | 204,964 | 56,993 | 2.57 | 9.25 |
| 768 | decaps | 644,420 | 236,201 | 141,976 | 2.73 | 4.54 |
| 1024 | keygen | 749,101 | 281,664 | 288,756 | 2.66 | 2.59 |
| 1024 | encaps | 803,207 | 276,782 | 74,342 | 2.90 | 10.80 |
| 1024 | decaps | 955,082 | 319,082 | 188,257 | 2.99 | 5.07 |

## Repository harness (`cargo bench --bench pqc_speed -- ml-kem`)

Base and new binaries alternated three times (`raw/pqc_speed_*.txt`); these
rows include OS randomness and 25-sample medians, so they are noisier than
the tables above and are kept as the cross-check, not the headline.
