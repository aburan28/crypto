# Stage 20 Phase B: terminal balanced-PDP results

The blinded production panel completed and was scored only after the terminal
run inventory was sealed. It contains 160 public synthetic instances and 480
backend outcomes across four balanced cells. Each cell has 20 decomposable and
20 nondecomposable targets. The source-system clustering check found 160
clusters of size one, so no target pair shares a source system.

Post-run truth scoring classified **261 outcomes as correct and 219 as
inconclusive**. There were **zero false positives and zero false negatives**.
Timeouts and conflict-capped `Unknown` outcomes remain inconclusive; they are
not counted as negative results.

The verified additive terminal evidence bundle under
`stage-20-phase-b-terminal-evidence-successor-04-20260910/` closes three
custody omissions in the original score format without altering that frozen
score. It binds the raw run-seal file hash, the scorer source and metered
command/output receipt, and the final WDSat build receipt alongside the Rust
and CryptoMiniSat builds. It also retains the public blind input, both frozen
protocols, all 21 campaign-level accounting receipts, and a self-hashed file
inventory. The report, summary, scorer, process meter, and bundle verifier are
archived inside the bundle, so later checkout changes do not invalidate the
measured source binding. The sealed oracle ledger itself is not copied.

Run `python3 scripts/verify_stage20_phase_b_terminal_evidence.py` to verify the
bundle, reconstruct its outcome totals and campaign accounting, and enforce the
unchanged non-SOTA claim boundary.

The Phase-B implementation and execution plan are bound to clean commit
`62be8cf6dfdaf1aec54e0a8f290ed55b822229a0`. Phase A was prepared from clean
commit `47235e51b74a6fa8f3c8dc85d68bf886e20a1e88`.

## Per-cell scored outcomes and charged backend processes

`Correct` is shown as total with `(true positive / true negative)`. Every row
has zero false positives and zero false negatives. Core time and summed process
wall include the main solver process and, for WDSat and CryptoMiniSat SAT
outcomes, the conditional external point-witness validator. The native-XOR
process already includes source regeneration, solving, source-model validation,
and exact lifted point-witness validation.

| Cell | Backend | Correct (TP/TN) | Inconclusive | Core-s | Summed process wall-s | Peak RSS MiB | Conflicts reported | Conflict sum | Median reported conflicts |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `n31-l5-m3-standard-a1-f0` | native-XOR | 40 (20/20) | 0 | 18.782658 | 20.050319 | 14.70 | 40/40 | 696,384 | 21,541.5 |
| `n31-l5-m3-standard-a1-f0` | WDSat | 40 (20/20) | 0 | 6.937482 | 8.070513 | 10.19 | 40/40 | 791,042 | 26,683 |
| `n31-l5-m3-standard-a1-f0` | CryptoMiniSat | 40 (20/20) | 0 | 48.107357 | 52.778182 | 13.58 | 40/40 | 3,460,100 | 104,966.5 |
| `n31-l5-m3-ggmp-a0-f0` | native-XOR | 2 (2/0) | 38 | 81.057445 | 89.968778 | 51.81 | 40/40 | 3,814,284 | 100,000 |
| `n31-l5-m3-ggmp-a0-f0` | WDSat | 1 (1/0) | 39 | 4,246.896465 | 4,701.108792 | 11.72 | 1/40 | 6,121,750 | 6,121,750 |
| `n31-l5-m3-ggmp-a0-f0` | CryptoMiniSat | 17 (17/0) | 23 | 2,967.797629 | 3,304.747181 | 285.77 | 17/40 | 16,996,950 | 839,886 |
| `n41-l5-m3-standard-a1-f0` | native-XOR | 40 (20/20) | 0 | 20.508182 | 23.843650 | 13.98 | 40/40 | 723,064 | 23,457 |
| `n41-l5-m3-standard-a1-f0` | WDSat | 40 (20/20) | 0 | 7.407406 | 8.888666 | 10.25 | 40/40 | 686,457 | 23,672 |
| `n41-l5-m3-standard-a1-f0` | CryptoMiniSat | 40 (20/20) | 0 | 16.399481 | 19.706110 | 10.83 | 40/40 | 1,264,702 | 36,696.5 |
| `n59-l9-m3-standard-a1-f0` | native-XOR | 0 (0/0) | 40 | 425.686974 | 487.320580 | 163.38 | 40/40 | 4,000,000 | 100,000 |
| `n59-l9-m3-standard-a1-f0` | WDSat | 1 (1/0) | 39 | 4,237.965518 | 4,779.780599 | 13.91 | 1/40 | 3,307,792 | 3,307,792 |
| `n59-l9-m3-standard-a1-f0` | CryptoMiniSat | 0 (0/0) | 40 | 4,358.053157 | 4,800.585750 | 360.48 | 0/40 | 0 | unavailable |

All three backends completely resolved the standard n=31 and n=41 cells. The
n=31 GGMP cell was substantially harder: CryptoMiniSat resolved 17 positive
instances, native-XOR resolved two, and WDSat resolved one; every negative
instance remained inconclusive. At n=59, native-XOR reached its conflict cap on
all 40 targets, both CryptoMiniSat classes reached their watchdogs, and WDSat
resolved one positive target. Those capped outcomes establish no SAT or UNSAT
answer.

Conflict medians above use only rows that reported a numeric conflict count.
No value is imputed for a timeout without a retained numeric count. WDSat thus
has one reported value in each of the GGMP n=31 and standard n=59 cells, while
CryptoMiniSat has no reported conflict value at n=59.

## Backend totals

| Backend | Correct (TP/TN) | Inconclusive | Core-s | Summed process wall-s | Peak RSS MiB | Conflicts reported | Conflict sum | Median reported conflicts |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| native-XOR | 82 (42/40) | 78 | 546.035259 | 621.183327 | 163.38 | 160/160 | 9,233,732 | 29,265.5 |
| WDSat | 82 (42/40) | 78 | 8,499.206871 | 9,497.848570 | 13.91 | 82/160 | 10,907,041 | 30,712 |
| CryptoMiniSat | 97 (57/40) | 63 | 7,390.357624 | 8,177.817224 | 360.48 | 97/160 | 21,721,752 | 51,280 |

The three backend totals sum to 16,435.599754 core-seconds and
18,296.849121 summed child-process wall-seconds. Including exports, isolated
source verification, and the other non-backend children raises the terminal
run-summary totals to 16,443.830358 core-seconds and 18,306.751722 summed
child-process wall-seconds. The separately retained whole-driver receipt is
inclusive: 16,494.118381 core-seconds, 18,373.393590 seconds elapsed wall, and
377,995,264 bytes maximum child high-water RSS. It encloses the child work and
must not be added to those child sums.

## Builds and successful terminal-path accounting

Final build receipts are bound into the execution plan. Their internal process
sums remain useful stage detail, but the campaign total below uses each
inclusive outer builder receipt once.

| Successful-path component | Core-s | Receipt wall-s | Peak RSS MiB |
|---|---:|---:|---:|
| Phase-A preparation | 1,072.186495 | 1,109.610274 | 2,081.02 |
| Sanitized Phase-B staging | 0.150859 | 0.537637 | 24.12 |
| Final Rust exporter/backend outer build | 145.970568 | 98.087391 | 894.47 |
| Final CryptoMiniSat outer build | 149.271767 | 97.291123 | 1,004.75 |
| Final WDSat capacity-successor outer build | 1.809385 | 2.250141 | 55.25 |
| Terminal Phase-B run outer receipt | 16,494.118381 | 18,373.393590 | 360.48 |
| Post-run scoring | 0.694939 | 1.282407 | 46.86 |
| **Successful terminal path** | **17,864.202394** | **19,682.452562** | **2,081.02 maximum** |

The wall total is an arithmetic sum of seven separate receipt intervals, not
elapsed calendar time. The RSS total is the maximum individual high-water mark,
not a sum. The Rust and CryptoMiniSat builders used two jobs; their recorded RSS
is the largest child high-water value and does not measure simultaneous
aggregate build-tree memory. Measured single-core elapsed time is unavailable.
Every legacy `single_core_seconds` field is only an alias for total user plus
system CPU.

## Retained failed and superseded work

The broader campaign accounting retains 21 non-nested top-level receipts:
15 returned zero and six returned nonzero. Adding the selected terminal path,
failed attempts, superseded clean builds, planning preflights, and the separate
dependency-fetch receipts gives **18,546.173826 core-seconds** and
**20,239.841542 summed receipt wall-seconds**. The maximum individual RSS
remains the 2,182,103,040-byte Phase-A preparation receipt.

| Non-terminal receipt | Disposition | RC | Core-s | Wall-s |
|---|---|---:|---:|---:|
| Initial plan preflight | failed aggregate-resource validation | 1 | 0.303718 | 0.496837 |
| Plan preflight successor 01 | superseded | 0 | 0.350476 | 0.595054 |
| Plan preflight successor 02 | selected supporting preflight | 0 | 0.343995 | 0.452549 |
| Initial Rust build | failed vendor step | 2 | 1.112891 | 2.357245 |
| Rust build successor 01 | superseded clean build | 0 | 199.908939 | 201.413121 |
| Rust build successor 02 | superseded clean build | 0 | 141.788862 | 96.281341 |
| Initial Rust dependency fetch | failed filesystem permission | 101 | 0.110110 | 0.534628 |
| Rust dependency fetch successor 01 | supporting prerequisite | 0 | 0.148980 | 0.376952 |
| Initial CryptoMiniSat build | failed binary-copy step | 2 | 149.689779 | 107.512258 |
| CryptoMiniSat build successor 01 | superseded clean build | 0 | 154.340930 | 104.288698 |
| Initial WDSat outer build | failed builder validation runtime | 1 | 2.408927 | 3.174999 |
| WDSat build successor 01 | superseded clean build | 0 | 1.845063 | 2.431117 |
| WDSat build successor 02 | superseded by capacity successor | 0 | 1.603593 | 2.345087 |
| Initial Phase-B run | retained capacity failure; no solver started | 2 | 28.015169 | 35.129094 |

The initial run exported and independently source-checked all 160 blind
instances, then stopped before any solver process. Its frozen sizing record found
that WDSat required `max_buffer_size=32804`, while the bound build supplied
`32264`. The capacity successor changed that limit through an additive protocol
and fresh WDSat build. The failed run contributes operational cost only and no
SAT/UNSAT result.

## Custody and identities

- Phase-A protocol: `69544ee3a0b0250dd720e1b498e8f03cdf7d4427ab11bfb551325cf29b86fd53`
- Phase-A seal: `01b7fbab8bccbb45c7e3e6a46640ff5997cefee303416afd9671da7cc2e13632`
- Blind bundle: `b45d0849e5ff4d3126e8dec1148b7f4da54e8a92da7b9b67d6be58700b7e9d77`
- Phase-B capacity-successor protocol: `5de2d43527286054ad348b182145d07a988b1236a141377d01b0e7f6f5a38ddb`
- Terminal run inventory: `bb4377054700edf3049a7590b9b9fdd8fe6a0c57b6a216a7928dd211dc74a12c`
- Terminal `run-seal.json` file: `e528f9bbaa6b7554224379cb32f7b1ed5be0a3cb8bde15e56dee7d1a5f75b538`
- Terminal run summary: `198a3ec9a4569f1a73a51ea92fa3d525475af038223408a3b1ca2de780459003`
- Scored result: `cc6a999e4438a01b44c187039f1fa7ff684d2c8f38192d7c4032159d62bf806e`
- `score-seal.json` file: `675f951516ecd9a81ac11f669dc8bf56a321d5da0a59b60cb04b784a7f0e7cab`
- Rust build receipt: `de34841bbf582143a41ae01d50b6fc5e029626e7c09cd254b895a5c1d652cb9a`
- CryptoMiniSat build receipt: `c2a9be07fb510e378b935433c8c091dd4f515e2f05e24154f8ef00e83fc191f4`
- WDSat capacity-build receipt: `7239950c441ca46554bbabdbd3c9e4768adfa6b2c784bd1b6696446b22054fa7`

The selected executable hashes are
`48c63078f2da6e594fef0038c6ecc13c34beb569f35660b6c39532123a2358e2`
for the exporter,
`14ec9a3d976371b6b100fa32a4ab916f50ddc88aad18fb6c6df8e81de925f96d`
for the isolated backend,
`6c509f09622f103d8a3ad90afc151e1c4031275c052d7f8af481465b4f27f2af`
for CryptoMiniSat, and
`1de1e3328af3f00963cd5339c07cc4ab18efd55d1749c953ecaafcefc54fd185`
for WDSat.

## Evidence boundary and remaining gates

The frozen score retains `full_cost_gate_passed: false` and
`independent_external_reproduction_satisfied: false`. This additive report now
lists Phase-A preparation and the final build-driver outer receipts, but it does
not mutate that score or close the remaining accounting and scientific gaps.

- No receipt measures true single-core elapsed time. Parallel build receipts do
  not measure simultaneous aggregate process-tree memory, and prior source,
  toolchain, system dependency, and licensing costs remain outside the bound
  build totals.
- The balanced 50/50 allocation measures conditional solver behavior. It does
  not estimate natural decomposability prevalence or relation yield.
- Phase B performs no relation collection, rank accumulation, modular linear
  algebra, unknown-scalar recovery, or same-target automorphism-optimized
  Pollard-rho comparison.
- Magma F4 was not executed in Phase B. There is no matched process-scoped
  Magma CPU/RSS record or Magma-linked point-witness evidence.
- There is no unaffiliated external reproduction or completed source-pinned
  novelty review.

The supported conclusion is therefore narrow: on this frozen balanced public
toy panel, every terminal answer from native-XOR, WDSat, and CryptoMiniSat was
correct, while all other outcomes remained inconclusive. This is not a natural
yield result, an end-to-end index-calculus result, a rho crossover, a novelty
finding, a key-recovery result, or a Koblitz index-calculus state-of-the-art
claim.

The compact machine-readable companion is
`stage-20-phase-b-result-summary-20260910.json`.
