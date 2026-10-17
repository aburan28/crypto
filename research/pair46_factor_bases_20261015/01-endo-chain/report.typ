#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#set heading(numbering: "1.")
= 01-endo-chain: factor bases closed under the degree-5 CM endomorphism (pair 46)
_2026-10-15. Pair-46 factor-base campaign, item 01. Protocol `PROTOCOL.md` (sha256 0ce98a09…dfcf2), frozen before data. Instrument: standalone Rust binary `src/main.rs`, 1 core, L0 timing. All 1 600 phase-1 solves and 800 post-hoc solves were verified by d·P = Q._

*Verdict: negative, with an exact reason (E3).* The best endo-chain arm, a single chain from P, costs
$S = 1.47 thin [1.28, 1.69]$ group-op equivalents per $sqrt(N)$ on the 22-bit toy. Strided ±BSGS
costs $0.955 thin [0.86, 1.05]$ on the same targets (ratio $1.54 thin [1.30, 1.82]$), and the rho-with-negation walk
takes $0.88 sqrt(N)$ steps. Chains do beat the abscissa (small-x) base, by a factor of about $sqrt(k)$ at fixed $|"FB"|$.
That is a constant, and it ends at the single-chain search, which stays above BSGS.

== Setting
Pair 46: $p = 5947698030651343$, $N = 5947697996161643$, $p = 5t^2+19t+19$, $t^2-4p = -19(t+2)^2$.
The endomorphism $w = (1+sqrt(-19))/2$ has degree 5 and eigenvalue $lambda = 2973849084305079$. Its x-map
$R = "Num"_5 / "Den"_4$ was interpolated and verified on 1000 points, and $"ord"(lambda) = (N-1)/7 = 849671142308806$.
Toy curves come from the same family: the crater curve $j = -96^3$ with prime order at 16, 18, 20 and 22 bits.

== Exact statements
- *E2.* There is one unknown per chain: $log psi^i H = lambda^i log H$. Two hits on the same chain give a linear equation in $d$.
- *E3.* A known-head chain with an x-only table is BSGS with table entries costing $c_E$ each and unstrided giant steps. Its optimum is $sqrt(2 c_E N)$ against $sqrt(N)$ for strided ±BSGS. More chains with unknown heads add birthday waiting, so one chain is optimal.
- *E4.* $"Aut" = {plus.minus 1}$, so glv-orbit reduces to the ± fold. Folding the walk by $psi$-orbits needs canonical representatives of orbits of length $8.5 dot 10^14$, which is infeasible.
- *E5.* Phase-1 artefact: with the fixed-$b$ walk $Q + a P$, the two hit positions of an element are symmetric about $-d$, so the abscissa arm costs $Theta(N)$ (measured $A/N = 0.236$). A post-hoc random-relation walk removes this.

== psi versus one addition (measured, ns per operation, this host)
Pair 46: addition 462.4 ns, psi x-evaluation 454.2 ns, so $c_E = 0.982$. Legendre symbol: 497.2 ns. Toys: $c_E$ = 1.06–1.14. Both operations are dominated by one inversion. In field operations, an affine addition is 1I+2M+1S and psi is 1I+10M. With batched inversion psi would cost about 2.2 additions (derived), so the measured $c_E approx 1$ favours endo-chain.

#figure(image("fig_cost.svg", width: 100%), caption: [Left: phase 1 (preregistered), mean $S$ with bootstrap 95% CI, $n = 40$ per arm and size (`results.json`). Right: post-hoc random-relation walk, $S$ against chain length $k$ at $F = ceil(sqrt(N))$ (`posthoc_v2_randomwalk.json`).])

== Phase 1 headline values (S, mean [95% CI], 22-bit toy, N = 2699659)
#table(columns: 2, [rho (set-up included; walk alone 1.25)], [2.02 [1.83, 2.24]], [rho-neg (restarts included; walk alone 0.88)], [2.93 [2.61, 3.26]],
[bsgs-pm], [0.955 [0.86, 1.05]], [abscissa-k1 (E5 artefact)], [390 [312, 463]], [endo-chain-k16], [8.32 [7.20, 9.59]],
[endo-chain-kF (unknown head)], [2.12 [1.90, 2.34]], [endo-chain-known], [1.47 [1.28, 1.69]], [k16-nofold (positive control)], [341 [272, 407]])
Endo-chain-known / BSGS by size: 1.49, 1.41, 1.63, 1.54, all with CI above 1.18. The E3 prediction is $sqrt(2 c_E) approx 1.46$–$1.51$.
The fold's positive control (nofold / k16) ranges from 14 to 41, so the instrument sees the fold.
Post-hoc ratio to abscissa at 22 bits: k16 0.29 [0.23, 0.36]; kF 0.096 [0.083, 0.114].

#figure(image("fig_m2_rate.svg", width: 80%), caption: [Phase 2 (preregistered): probability that a random point is $plus.minus P_i plus.minus P_j$, Wilson 95%, $T = 3000$, 20-bit toy, $|"FB"| = 201$. The endo-chain sets fall inside the random-set null. The arithmetic-progression positive control is about 100 times lower.])

== Limits
Toy sizes only (16–22 bits; pair 46 was used for E1 and the timing). The host is shared, so timings are not isolated. Decomposition was tested for $m = 2$ only. This is not an ecbench session. Nothing was cut from the preregistered design, and the post-hoc arms are labelled. Reproduce with the commands in `RESULTS.md`.
