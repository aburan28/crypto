#set page(paper: "a4", margin: 1.6cm)
#set text(size: 9pt)
#set table(stroke: 0.4pt)
#show raw: set text(size: 8pt)
#show figure.caption: it => it.body
#show table: set text(size: 7.8pt)

= The root-finding decomposition oracle inside ecbench's index calculus: binary and prime ladders against rho

2026-10-13. Studies `research/pdp_root_finder_ecbench_20261013/` (binary ladder, its claim round in `claim_round/`) and `research/pkm_tower_ecbench_20261013/` (PKM's prime-field Kummer tower). Oracles: `ic_framework::rootfind` (`rootfind-mitm`, `rootfind-pairs`). Factor base for prime curves: `ic_framework::prime_tower` (`prime-kummer`). Protocols registered before the runs; one amendment (the oracle's accounting, before the binary ladder restarted). *Class: engineering.*

== What was measured

Index calculus end to end — factor base, pair table, relations, linear algebra, verification — in ecbench's one unit, $S$ = charged group-addition equivalents $\/ sqrt(r)$, on one-target workloads, every run verified by the harness and every run replayed by the audit.

- *Binary ladder.* $E_0: y^2 + x y = x^3 + 1$ over $FF_(2^n)$, $n = 13, 19, 23, 31, 37$, standard subspace base of dimension $d = 5, 7, 8, 10, 12$; arms: strong signed-Frobenius rho (reference), the inherited-F4 Gröbner oracle (`pdp3-koblitz`, $n <= 23$), the harness's `mitm:m=3` table oracle, `rootfind-mitm:m=3`, `rootfind-pairs:m=3` ($n <= 31$), an A/A rho control. Eight public targets, five rounds ($n <= 19$) or three.
- *Prime ladder.* $y^2 = x^3 + a x + b$ over $FF_p$, $p = c dot 2^12 + 1$ of 20–28 bits, PKM's base $V = mu_(2^t)$, $t = 8 dots 11$; arms: rho with the negation map (reference), `mitm`, `rootfind-mitm`, `rootfind-pairs` ($t <= 9$), the same oracle on a plain abscissa base of equal size (control), A/A. Eight targets, three rounds.
- *Claim round.* `rootfind-mitm` on a Frobenius-orbit base (`koblitz-orbit`) at $n = 23, 31$, so that a `vs_rho` claim can be built; one claim per workload.

All 1,248 binary-ladder runs, 896 prime-ladder runs and 256 claim-round runs verified $[d]G = Q$; every replay audit is `OK` with every replay identical.

== Binary ladder

#figure(image("figures/s_by_size.svg", width: 100%), caption: [Every measured run's $S$ by arm and size, binary ladder, with the floor $sqrt(pi\/4n)$. The F4 arm's $S$ is a lower bound: its Gröbner word XORs are counted and not charged.])

#table(columns: 7, align: (right, right, right, right, right, right, right),
  [$n$], [`rootfind` ÷ strong rho, $S$], [÷ strong rho, wall], [÷ F4, $S$ (F4 a lower bound)], [÷ F4, wall], [÷ `mitm`, $S$], [÷ `mitm`, wall],
  [13], [2.85 [2.64, 3.04]], [4.18], [2.56], [*0.0029* (admitted)], [1.263], [0.953],
  [19], [15.0 [14.3, 15.7]], [11.6], [9.59], [0.0004 (descriptive)], [1.186], [0.697],
  [23], [60.9 [56.1, 66.0]], [34.5], [29.9], [0.0001 (descriptive)], [1.047], [0.427],
  [31], [3,053 [2,505, 3,596]], [2.33], [—], [—], [1.003], [0.278],
  [37], [25,116 [20,793, 30,546]], [descriptive], [—], [—], [1.0008], [0.312 (descriptive)],
)

Paired totals ratios from `ecbench compare` with 95 % cluster-bootstrap intervals; wall ratios are the harness's admitted medians where both runs held L2. Against the F4 oracle, the whole solve is 350× faster in wall time at $n = 13$ (admitted) and about 2,500× and 10,000× at $n = 19, 23$ (descriptive; the F4 runs are too long to hold L2). Against rho, index calculus costs more on every size and the gap grows; the jump at $n = 31$ is mostly the curve's cofactor (1,492): index calculus pays for the whole group, rho only for the prime subgroup. The relations phase is 71–94 % of $S$, the pair table 6–29 %, factor base and linear algebra under 2 %.

#figure(image("figures/phase_split.svg", width: 100%), caption: [Mean $S$ by phase for the binary ladder's index-calculus arms (verified runs).])

== Prime ladder (PKM's Kummer tower)

#figure(image("../pkm_tower_ecbench_20261013/figures/s_by_size.svg", width: 100%), caption: [Every measured run's $S$ by arm and size, prime ladder, with the floor $sqrt(pi\/4)$. Every arm charges every group operation.])

#table(columns: 6, align: (right, right, right, right, right, right),
  [bits of $p$], [$log_2 r$], [`rootfind` ÷ rho, $S$], [÷ rho, wall], [tower ÷ plain base, $S$], [`rootfind` ÷ `mitm`, $S$ / wall],
  [20], [18.85], [34.0 [27.9, 42.6]], [15.5], [1.006], [1.21 / 1.31],
  [22], [21.19], [41.9 [34.3, 53.4]], [26.4], [1.001], [1.14 / 1.13],
  [24], [23.38], [29.4 [23.7, 37.3]], [22.8], [1.012], [1.08 / 1.17],
  [26], [23.28], [96.6 [77.8, 124.1]], [83.4], [1.005], [1.06 / 1.06],
  [28], [26.20], [140.6 [113.0, 182.3]], [178.6], [0.999], [1.03 / 1.09],
)

Index calculus with PKM's base costs 29–141× rho with the negation map, growing by $0.26 plus.minus 0.28$ bits per bit of $r$ (least squares over the five sizes), consistent with the $p^(1\/6)$ the floors predict and too wide here to separate from zero. The tower base costs the same as a plain base of equal size within 1.2 %: its structure helps an algebraic solver, not a table. On prime curves the new oracle is 3–31 % slower than `mitm` in wall time: the prime group has no batched addition, so the block of 32 abscissae only adds overshoot.

== Claims

On the binary ladder's subspace base every `vs_rho` claim build was refused: the IC1 identity names a base by its sign-and-Frobenius orbits and that base has more columns than orbits. On a Frobenius-orbit base, rerun with the binary that includes main's claim-builder fix, all 16 claims build with `operation_accounting` and `record_class` and are checked. Every claim shows index calculus slower than the strong rho on its target (online speedup 0.132–0.232 at $n = 23$, 0.005–0.017 at $n = 31$). Every claim fails the checker on independent validation (no replay receipt from a second hardware class is available in this container); 11 are classed exploratory because a run of the pair held L1 instead of L2.

== Decision

The root-finder replaces the Gröbner oracle inside the index-calculus pipeline at thousands of times less wall time on binary curves, matches the harness's own table oracle in charged work, and leaves the method's growth unchanged: on binary and prime curves alike, index calculus with an enumerating oracle stays above rho and the gap widens with the curve, as the pair-table floor ($p^(2\/3)$ at $m = 3$) against rho's $p^(1\/2)$ says it must. Moving that exponent needs a decomposition test whose per-trial cost does not grow with the factor base; the scope note `research/ic_speed_and_rho_scope_20261013/` lists the candidates and their falsifiers.

Environment: Claude Code cloud container, Intel Xeon \@ 2.10 GHz, 4 vCPU, L3 260 MiB, 15 GiB RAM, kernel 6.18.44, rustc 1.97.0; `ecbench run --cpus auto` (a reserved core, other threads evicted), one thread. Conductor's control plane was unreachable in this container: no task card. Canonical graphs checked and left unchanged: the IC scoreboard, `progress-timeline.json`, the leaderboard and `docs/bounds/FRONTIER.md` (no result here changes a frontier: every ratio to rho is above 1).
