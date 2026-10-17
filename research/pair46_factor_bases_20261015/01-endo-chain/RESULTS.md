# 01-endo-chain — results

Campaign: pair-46 factor bases (2026-10-15), item 01. Protocol: `PROTOCOL.md`
(sha256 `0ce98a09…dfcf2`, frozen before data; recorded in `results.json`).
Instrument: standalone Rust binary `src/main.rs` (this directory), 1 core (`taskset -c 0`).
Total compute: phase 0–2 run 22.6 s, post-hoc runs 4–5 s each. Nothing was cut.

**Verdict: negative (with an exact reason).** The best endo-chain arm costs
S = 1.47 [1.28, 1.69] group-op equivalents / sqrt(N) on the 22-bit toy, against
+-BSGS 0.955 [0.86, 1.05] on the same targets (ratio 1.54 [1.30, 1.82]) and the rho-with-negation walk
0.882 sqrt(N) (measured walk steps; theory 0.886). Chains do beat the abscissa (small-x) base,
by a factor that depends on chain length k but not on the field, so that gain is a constant.

## Exact statements (verified computations / identities)

| id | statement | evidence |
|:--|:--|:--|
| E1 | On pair 46 the x-map of the endomorphism with eigenvalue lambda = 2973849084305079 is R = Num/Den, Num degree 5, Den monic degree 4 (coefficients in `results.json` → `pair46.xmap_num/den`). lambda^2 - lambda + 5 = 0 mod N; N P = O; p = 5t^2+19t+19, N = p+1-t | interpolated from 10 points, then checked x(lambda P) = R(x(P)) on 1000 random points (rational maps of degree ≤ 5 that agree at > 10 points are equal) |
| E1' | ord(lambda) = 849671142308806 = (N-1)/7 (matches the task header) | exact order from the factorisation of N-1 |
| E2 | One unknown per chain: log psi^i(H) = lambda^i log H. Two hits (a1+b1 d) = s1 lambda^{i1} u and (a2+b2 d) = s2 lambda^{i2} u on the same chain eliminate u and give a linear equation in d | algebra; every one of the 1 600 phase-1 solves and 800 post-hoc solves was checked by d P = Q |
| E3 | A single chain from a known head (P itself) with an x-only table is a baby-step giant-step search. Each table entry costs c_E group-op equivalents, and the giant steps are random rather than strided. So the optimum is sqrt(2 c_E N) against sqrt(N) for strided +-BSGS, a factor sqrt(2 c_E) ≥ 1.41 whenever c_E ≥ 1. More chains with unknown heads only add birthday waiting among C chains (cost grows like sqrt(C)), so k = F (one chain) is optimal and the method cannot fall below BSGS | cost identity; matches the measurements (known-chain/BSGS = 1.41–1.63) |
| E4 | glv-orbit needs a non-trivial automorphism group. j = -96^3 ≠ 0, 1728 gives Aut = {±1}, so glv-orbit reduces to the ± fold, which the abscissa arm already has. A psi-orbit fold of the walk (the glv-orbit trick applied to psi) would need canonical representatives of orbits of length ord(lambda) = 8.5e14 on pair 46. That is not computable per step | exact |
| E5 | Phase-1 artefact: the fixed-b walk Q + aP can only hit an unknown element u at a = u - d and a = -u - d. These two positions are symmetric about -d (mod N), so a same-element pair needs the walk to reach -d or -d+N/2, which takes Θ(N) steps. Observed: abscissa-k1 spent A/N = 0.236 at 22 bits. Distinct exponents break the symmetry for chains, so the phase-1 abscissa-vs-chain ratio is inflated. The post-hoc random-relation instrument (below) removes this | exact argument + measurement |

## Cost of psi versus one addition (pair 46 and toys) — measured

| curve | ns / affine add | ns / psi x-eval | c_E = psi/add | ns / Legendre (c_L) | field ops |
|:--|--:|--:|--:|--:|:--|
| pair 46 (53-bit) | 462.4 | 454.2 | **0.982** | 497.2 (1.075) | add 1I+2M+1S; psi 1I+10M |
| toy 16 / 18 / 20 / 22 | – | – | 1.069 / 1.058 / 1.143 / 1.063 | (0.75–0.86) | same |

Both operations are dominated by the one inversion (u128 extended Euclid, about 450 ns), so c_E ≈ 1
here. Under a batched-inversion model (I ≈ 3M) an addition costs about 6M and psi about 13M, so c_E ≈ 2.2.
The measured c_E ≈ 1 therefore favours the endo-chain arms (derived, not measured). Phase 1
also reports S with c_E = 1 (`S_cE1_*`); the conclusions do not change.

## Phase 1 (preregistered): full DLP solves, 40 targets per arm per size, all verified 40/40

S = (A + c_E E + c_L L)/sqrt(N), mean with bootstrap 95% CI (measured). F = ceil(sqrt N).

| arm | 16-bit | 18-bit | 20-bit | 22-bit |
|:--|:--|:--|:--|:--|
| rho (r=20, full storage) | 5.09 [4.90,5.28] | 3.67 [3.50,3.86] | 2.51 [2.32,2.71] | 2.02 [1.83,2.24] |
| rho-neg | 5.56 [5.30,5.83] | 4.27 [3.91,4.67] | 3.44 [3.07,3.85] | 2.93 [2.61,3.26] |
| bsgs-pm | 1.17 [1.08,1.26] | 1.17 [1.08,1.25] | 1.04 [0.95,1.13] | **0.955 [0.86,1.05]** |
| abscissa-k1 (small-x) | 54.7 [44.7,64.5] | 107 [88,124] | 205 [171,239] | 390 [312,463] |
| endo-chain-k4 | 6.36 [5.60,7.09] | 9.53 [8.34,10.8] | 11.2 [9.8,12.5] | 17.2 [14.0,21.0] |
| endo-chain-k16 | 4.37 [3.87,4.93] | 4.88 [4.27,5.52] | 5.14 [4.44,5.88] | 8.32 [7.20,9.59] |
| endo-chain-k64 | 2.75 [2.49,3.02] | 3.45 [2.96,3.97] | 3.70 [3.21,4.23] | 4.61 [4.11,5.11] |
| endo-chain-kF (1 chain, unknown head) | 2.20 [2.08,2.34] | 2.43 [2.18,2.72] | 2.35 [2.14,2.58] | 2.12 [1.90,2.34] |
| endo-chain-known (1 chain from P, k*) | 1.75 [1.51,2.02] | 1.64 [1.40,1.90] | 1.69 [1.50,1.91] | **1.47 [1.28,1.69]** |
| endo-chain-k16-nofold (positive control) | 60.9 [50.5,70.8] | 119 [99,139] | 191 [157,225] | 341 [272,407] |

Rho rows include fixed set-up costs: 20 random multipliers (40 scalar multiplications), the start
point and, for rho-neg, a full restart after each fruitless cycle, all charged. The walk steps alone are
1.25 / 1.27 / 1.10 / 1.25 sqrt(N) for rho (theory 1.253) and
0.89 / 0.88 / 0.86 / 0.88 sqrt(N) for rho-neg (theory 0.886) (`walk_steps_over_sqrtN`). The
set-up is a constant number of group operations, so rho's S tends to its walk value as N grows.
Read the rho rows as an upper bound on the reference.

Derived ratios (bootstrap ratio of means, `derived_ratios.json`):
- endo-chain-known / bsgs-pm: 1.49 [1.26,1.77], 1.41 [1.19,1.68], 1.63 [1.41,1.88], 1.54 [1.30,1.82]. Predicted by E3: sqrt(2 c_E) = 1.46–1.51.
- endo-chain-kF / bsgs-pm: 1.88, 2.09, 2.26, 2.22 (all CIs above 1.7).
- Positive control for the fold (nofold / k16): 14.0 [11.1,17.0], 24.4, 37.2, 40.9 [31.7,52.4]. The instrument sees the fold.
- k16 / abscissa-k1 (phase 1): 0.080 down to 0.021. This is inflated by artefact E5.

## Post-hoc supplementary (NOT preregistered): random-relation walk

The walk R ← R + M_j uses r = 20 multipliers c_j P + e_j Q, with j drawn from an RNG and not from x.
Hits are decoded with general (a, b). F = ceil(sqrt N) is fixed and only k varies
(`posthoc_v2_randomwalk.json`):

| k | 16-bit | 18-bit | 20-bit | 22-bit |
|:--|:--|:--|:--|:--|
| 1 (abscissa) | 14.1 [12.7,15.7] | 17.4 [15.0,19.9] | 21.4 [18.7,24.2] | 28.9 [24.8,32.6] |
| 4 | 9.8 | 12.1 | 12.0 | 12.0 [10.6,13.6] |
| 16 | 7.5 | 7.7 | 7.2 | 8.3 [7.0,9.7] |
| 64 | 6.5 | 5.5 | 4.8 | 5.5 [4.8,6.1] |
| F | 6.1 | 4.4 | 3.5 | 2.77 [2.57,2.98] |

Ratio to abscissa at 22 bits: k4 0.42 [0.35,0.50], k16 0.29 [0.23,0.36], k64 0.19 [0.16,0.23],
kF 0.096 [0.083,0.114]. These fit the birthday-over-chains heuristic: at fixed F the gain is
about sqrt(k), so the fold is a constant set by k. With k = F the arm is the single-chain search of E3, and that
stays above BSGS. A first post-hoc variant
(`posthoc_v1_xdetermined.*`) chose the multiplier from x. That turns the factor base into a
distinguished-point set, so the arm becomes van Oorschot–Wiener rho
(S ≈ 6.9–12 for every k). It is kept as a record but not used.

## Phase 2 (preregistered): m = 2 decomposition structure, 20-bit toy, |FB| = 201, T = 3000

| factor base | rate P(R = ±Pi ± Pj) | Wilson 95% |
|:--|--:|:--|
| random points (null) | 0.1167 | [0.106, 0.129] |
| small-x | 0.1203 | [0.109, 0.133] |
| endo-chain k=16 | 0.1213 | [0.110, 0.134] |
| endo-chain single chain | 0.1170 | [0.106, 0.129] |
| arithmetic progression (positive control) | 0.0010 | [0.0003, 0.0029] |

The random-set heuristic 1 - exp(-2F²/N) gives 0.118. The positive control separates from the null
by about 100x, so the detector works. Both endo-chain sets fall inside the null interval: under the
psi-closure, sums of two factor-base elements behave like those of a random set of the same size.
No membership or decomposition structure was detected (comparative claim, one 20-bit curve, one FB size).

## Scope and limits
- Toys use t = smallest admissible with b-bit p in the pair-46 family p = 5t²+19t+19 (crater curve
  j = -96³ or its twist, prime order). ord(lambda) = (N-1)/18 and (N-1)/78 on the 16- and 18-bit toys,
  and (N-1)/2 on the 20- and 22-bit toys. On pair 46 it is (N-1)/7.
- Unit: affine group operations plus psi and Legendre operations converted by measured ns ratios on this host,
  which runs another experiment at the same time (one core, not isolated). Timings are L0 at best.
- Decomposition with m ≥ 3 and Semaev-style membership were not tested. E3 covers table-based use
  of chains. The phase-2 null result covers m = 2 only.
- Not an ecbench session; the instrument is a standalone binary in this directory.

## Files
`PROTOCOL.md`, `src/main.rs`, `Cargo.toml`, `run.log`, `results.json`, `posthoc_v2_randomwalk.{json,log}`,
`posthoc_v1_xdetermined.{json,log}`, `derived_ratios.json`, `ratios.log`, `plot.py` (layout only),
`fig_cost.svg/png`, `fig_m2_rate.svg/png`, `report.typ`, `report.pdf`.

## Reproduce
```sh
cd research/pair46_factor_bases_20261015/01-endo-chain
cargo build --release --offline
taskset -c 0 ./target/release/endo_chain 40 1800          # phases 0-2 -> results.json
taskset -c 0 ./target/release/endo_chain posthoc          # -> posthoc_v2_randomwalk.json
taskset -c 0 ./target/release/endo_chain posthoc xdet     # -> posthoc_v1_xdetermined.json
./target/release/endo_chain ratios                        # -> derived_ratios.json
python3 -m venv /tmp/rv && /tmp/rv/bin/pip install -q typst matplotlib
/tmp/rv/bin/python plot.py && /tmp/rv/bin/python -c "import typst; typst.compile('report.typ', output='report.pdf')"
```
The counts and seeds are deterministic. c_E and c_L come from wall-clock timing, so S values move
slightly between runs (for example, the 16-bit post-hoc abscissa S was 11.976 in a first run and 11.971 on rerun).
