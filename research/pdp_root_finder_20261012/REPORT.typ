#set page(paper: "a4", margin: 1.6cm)
#set text(size: 9pt)
#set table(stroke: 0.4pt)
#show raw: set text(size: 8pt)
#show figure.caption: it => it.body
#show table: set text(size: 7.6pt)

= Root-finding decomposition oracles for the point decomposition problem, measured on both tracks

2026-10-12. Study directory `research/pdp_root_finder_20261012/` (README: the full results document; `PROTOCOL.md`: the pre-registration; `reports/run-20261012T2135Z/`: the run report of the timing stage). Requirements and gates: `research/pdp_root_finder_requirements_20261012/` (PR \#1734). Driver: `examples/pdp_rootfind/`. *Class: engineering* — every arm's counted work sits on a generic floor (pair table or enumeration); the constants, the exactness instrumentation and the measured distance to F4 are the result; no exponent moved.

== What was built and measured

Two mechanisms meet the structural requirements S1–S8 without a Macaulay step: the memoryless *pairs-and-solve* by the group law (binary `rf-pairs`, prime `pairs3`; `≈ 10` multiplications a pair, $|F|^2\/2$ pairs a target) and the *2-sum table* (binary `rf-mitm`, prime `mitm3`, `mitm4`; keys $x(P_i plus.minus P_j)$, $5|F| - 3$ multiplications, one batched inversion and $2|F|$ probes a target, $|F|^2$ keys at 12 bytes). Probes read one 64-byte bucket of sixteen-bit fingerprints (no false negatives, about $16\/2^16$ false positives confirmed by binary search), all candidates of a target computed first and probed in one prefetched pass. References: PS-0's `decompose` (`baseline`) and PS-2's row kernel (`group-scalar`) on the binary rungs, `brute_rational` and the PKM checker's `algebraic_count` on the prime cells, the committed `f4_fp_tower` rows re-read (not re-run), the inherited F4 and the split-dreg F4 on a (23, 8) binary cell.

Cells: PS-0's frozen corpora $ell = 8 dots 12$ ($n = 3ell$, $y^2 + x y = x^3 + x^2 + 1$), fresh corpora at $ell = 13, 14$ and at $(n, ell) = (23, 8)$ from the same generator and seed rule, exhaustive digests; the Kummer-tower cells of the PKM rounds at $p_1 = 2013265921$ ($m = 2, 3, 4$), both targets of every committed row, and a fresh ladder $t = 4 dots 13$ ($m = 3$), $4 dots 12$ ($m = 2$), $4 dots 11$ ($m = 4$) with 16 uniform + 8 planted targets a rung.

== Exactness (G0)

Binary: verdict and root-set digests equal to the corpus's for both arms on every rung $ell = 8 dots 14$ (80 targets a rung at $ell <= 12$, 20 at $ell = 13, 14$), 0 mismatches; on (23, 8) the verdict digest is equal and 79 of 80 root sets compared equal, the 80th target being degenerate ($x_R in V$) and routed to the baseline decomposer by construction. Prime: 568 committed rows re-read; 707 / 707 arm root sets equal to the brute-force reference; 479 / 479 finished F4 verdicts equal to the algebraic count (18 unfinished F4 rows excluded); 0 of 486 rows with an algebraic solution and no rational decomposition; every planted target found on every ladder rung; 0 brute-force mismatches at $t <= 7$ ($m = 2, 3$) and $t <= 6$ ($m = 4$); verdicts of all arms equal at every $t$.

== Counted work: the one table

#figure(image("figures/binary_cost.svg", width: 100%), caption: [Figure 1. Binary rungs: wall time per refuted target on the reserved core (median of five isolated rounds; dots: rounds) and counted multiplications per target, with the fitted bits per rung. `rf-mitm` is the pair-table floor (1 bit per rung), `rf-pairs` the enumeration floor (2).])

Binary, $m = 3$, per refuted uniform target (medians; $T_"enum" = 2^(2 ell)\/2$ candidate pairs, pair floor $2^ell$ probes; wall columns from the five isolated rounds; `baseline` and `group-scalar` carry no counters):

#include "figures/table_binary.typ"

Prime, Kummer towers at $p_1$, full scan of 8 uniform targets (the cost of a refutation; $T_"enum"(m, t) = 2^((m-1)t)\/(m-1)!$, pair floor $2^(floor(m\/2) t)$):

#include "figures/table_prime.typ"

#figure(image("figures/prime_cost.svg", width: 100%), caption: [Figure 2. Prime ladder: counted multiplications per full scan against $t$ with the committed F4 rows (multiply-adds, finished runs) as context, and wall time per full scan on the reserved core.])

Amortised accounting (G5), from the recorded build counts: the 2-sum table costs about 5 multiplications a key, half of what the memoryless arm spends on one target, so it pays for itself at $K = 1$: binary $ell = 13$, $T_"prepare" = 41,032,579$ mul and $T_"target" = 20,257$, i.e. $T_"prepare"\/K + T_"target" = 4.1 dot 10^7$, $2.58 dot 10^6$, $1.81 dot 10^5$ at $K = 1, 16, 256$ against $8.21 dot 10^7$ for `rf-pairs`; prime $t = 13$: $4.17 dot 10^7$ to build against $8.34 dot 10^7$ for `pairs3`.

== Against F4 (G2)

#figure(image("figures/f4_vs_rootfind.svg", width: 100%), caption: [Figure 3. The committed `f4_fp_tower` rows (Kummer, $p_1$; multiply-adds, finished runs, medians of the rows of a cell) against the root-finders on the same cells and targets (field multiplications).])

Prime, same cells and targets as the committed rows (finished F4 runs only; a multiply-add and a multiplication compared one to one, which flatters F4):

#include "figures/table_f4.typ"

G2(a) asks for $<= 0.1 times$; the ratios are $10^5$–$10^12$ and grow with the cell. P4 held on every clause (`mitm3` 197 and `pairs3` 8,357 multiplications at $m = 3, N = 18$ against $8.99 dot 10^13$; `mitm4` 270 at $m = 4, N = 16$ against $2.71 dot 10^13$; every $m >= 3$, $N >= 15$ ratio above $10^9$). G2(b): the finished F4 rows grow by $3.214 plus.minus 0.063$ bits per unit of $t$ at $m = 2$ and $6.556 plus.minus 0.864$ at $m = 3$ (1.61 and 2.19 bits per tower variable, the figures the requirements document quoted), the root-finders by $1.05 plus.minus 0.09$, $1.02 plus.minus 0.02$ and $1.99 plus.minus 0.03$: disjoint intervals over $>= 4$ sizes on the prime track. Binary (23, 8): the inherited F4 decision oracle costs $4.10 dot 10^7$ word XORs and 231.5 ms a target (medians of 8, shared host) where `rf-mitm` costs 592 multiplications and 3.9 µs (isolated); the split-dreg F4 did not finish (killed, then an allocation failure under a 12 GB cap: censored). G2(c), the external tools on exported tower systems, was not run.

== Growth and the generic line (G3, G5)

#include "figures/table_slopes.typ"

Every root-finder's slope equals its floor within $plus.minus 1.96$ s.e. (P3, P5). The table arm's $S T^2$ gap to the generic line ($S = |F|^2$ keys, $T = ("mul" + "sqr")\/6$ group-operation equivalents, $r = 2^(n-1)$, $|"Aut"| = 2$) is positive and growing — that line prices a whole DLP, and a decomposition oracle alone cannot pass it; the row is the G5 accounting:

#include "figures/table_st2.typ"

*G3 fails; nothing here moves an exponent.* G6 (end-to-end budgets) follows: with the floors as the measured slopes, the 2-sum oracle at ECC2K-130, $m = 4$ ($|F| approx 2^33.9$) probes $approx 2^67.8$ keys a target against a budget of $2^35$; at $m = 5$, $2^81$ against $2^41$; at a 256-bit prime, $m = 4$, $t = 64$: $2^126$ against $2^59$. Extrapolation, not measurement.

== Wall time (reserved core, five rounds)

Medians of the uncontended rounds (`raw/bin_time.jsonl`, `raw/prime_time.jsonl`, `raw/isolated.jsonl`; 177 of 185 commands ran, 8 were refused by the settle check and leave fewer rounds on six $ell <= 12$ rows, 7 runs flagged contended are excluded). Binary, per refuted target: `rf-mitm` 3.1 µs at $ell = 8$, 57.6 µs at $ell = 12$, 156 µs at $ell = 13$, 482 µs at $ell = 14$; PS-0's `baseline` 9.3 ms, 3.07 s, 13.1 s, 55.5 s (3,007×, 53,255×, 84,095×, 115,160×); PS-2's `group-scalar` 276 µs, 48.1 ms, 191 ms, 757 ms (89×, 834×, 1,225×, 1,572×); `rf-pairs` within 0.91–1.18× of `group-scalar` on every rung (P2). The (23, 8) cell: 8.7 ms, 214 µs, 202 µs, 3.9 µs. `rf-mitm`'s wall grows by $1.196 plus.minus 0.123$ bits per rung against 1.000 in counted work: the per-probe cost rises from 12 ns ($ell = 8$, 0.2 MB table) to 29 ns ($ell = 14$, 805 MB, past the 260 MiB L3); the table build is recorded apart (8.2 s at $ell = 14$) and charged in §4. Prime, per full scan: `mitm3` 0.33 µs at $t = 4$ → 306.5 µs at $t = 13$, `pairs3` 2.3 µs → 574 ms (1,873× apart at $t = 13$), `subtract2` 0.23 → 61.7 µs ($t = 12$), `mitm4` 5.5 µs → 94.6 ms ($t = 11$); wall slopes 1.054 ± 0.089, 1.962 ± 0.031, 0.998 ± 0.070, 2.009 ± 0.071 bits per unit of $t$. A/A control (the baseline against itself at $ell = 10$, every round): 155–177 ms a target, at most 7.5 % apart within a pair. P1's wall clauses ($> 50 times$ at $ell = 8$, $> 1,000 times$ at $ell = 12$) held with 3,007× and 53,255×. Against the Gröbner references in wall: inherited F4 at (23, 8) 231.5 ms a target (shared host) against 3.9 µs, about $5.9 dot 10^4$ (unlike measurements, an order of magnitude); the committed $m = 3, N = 18$ F4 rows 10,968 s on their host against a 1.2 µs full scan here.

== Gates read, and status

#table(columns: 3,
  [gate], [result], [where],
  [G0 exactness], [pass on every arm and cell], [README §3],
  [G1 vs enumeration], [pass for the table arms (by memory, not algebra: ratio $7.8 dot 10^(-3) arrow 1.2 dot 10^(-4)$ binary, $0.11 arrow 2.4 dot 10^(-4)$ prime $m = 3$, halving a rung); fail (ratio ≈ 1) for the memoryless arms], [§4],
  [G2(a) vs F4, counted], [pass: prime $10^5$–$10^12$; binary: the Gröbner arm that finishes is about $6 dot 10^4 times$ slower in wall, the other does not finish], [§5],
  [G2(b) vs F4, growth], [pass on the prime track ($m = 2, 3$, disjoint intervals over ≥ 4 sizes); binary: F4's side is a two-size ratio, not a fit], [§5–6],
  [G2(c) vs msolve], [not run], [—],
  [G3 exponent], [fail: every slope is its floor, $S T^2$ gap positive], [§6],
  [G4 vs SAT], [not run], [—],
  [G5 accounting], [reported], [§4, §6],
  [G6 end-to-end], [fail (extrapolation at the floors)], [§8],
  [G7 class], [engineering on every row], [—],
)

Deviations from the protocol (README §2): the Gröbner reference cell is (23, 8), not (24, 8) (the F4 drivers refuse $n = 24$); the split-dreg reference is censored (memory); prime timing is the full scan of 8 uniform targets at every $t$ (at $t >= 12$ every uniform target decomposes); binary timing at $ell >= 13$ uses 4 refuted targets and one pass; the lookup structure changed before the timing stage with counted units unchanged and exactness repeated; slope intervals are $plus.minus 1.96$ s.e. over rungs (a bootstrap over targets has zero width for a table arm); the ladders and the committed rows read are supersets of the registered ones; R2, R4, R5, R8 not run.

On the request's wording ("the fastest algorithm you can make"): within the arms measured here the 2-sum table is the fastest on every cell, in counted work and in wall time, by the factors above; its exponent is the known pair-table floor, and nothing measured here bears on the end-to-end sizes. Claims stop there.

Environment: Claude Code cloud container, Intel Xeon \@ 2.10 GHz, 4 vCPU, L3 260 MiB, 15 GiB RAM, kernel 6.18.44, rustc 1.97.0, base commit `cd0286878`; timing through `tools/isolated_bench.py run --cpus 3` (CPU 3 reserved, other threads evicted, contended runs listed and not pooled), `RAYON_NUM_THREADS = 1`. Conductor's control plane was unreachable in this container (no Docker daemon or Postgres): no task card attached. Python lays out recorded values only; every table above is generated from `ANALYSIS.json`, written by the Rust `analyze` subcommand from the raw records. Canonical graphs checked and left unchanged: the IC scoreboard, `progress-timeline.json`, the leaderboard and `docs/bounds/FRONTIER.md` (no end-to-end $S$).
