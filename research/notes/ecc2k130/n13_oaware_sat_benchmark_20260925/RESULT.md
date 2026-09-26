# n13-m5 dense O-aware CNF: three cold SAT engines admitted on fixed toy corpus

The frozen [PR #785](https://github.com/aburan28/crypto/pull/785) panel
completed all 96 cold engine/target queries, with six of six tiny interface
smokes passing first. CryptoMiniSat 5.14.7, Kissat 4.0.4 and CaDiCaL 3.0.1
each returned **5 SAT and 27 UNSAT** over the same 32 exact Q+T targets.
Every SAT model selected one factor x per rotated slot; a separate
bit-serial/Fermat group implementation derived each complete rational
point fibre, tried at most 32 sign tuples, and found a signed tuple summing
to the exact target point. All 81 UNSAT statuses agree with the complete
#767 point oracle, but **no solver UNSAT proof certificate was independently
checked**. No target timed out, exceeded 2 GiB sampled child-tree RSS, or
returned an invalid model. Archive-only replay reconstructed every unit-
assumption CNF hash and all 96 raw statuses/witnesses, then passed.

| Engine | Exact SAT / oracle-consistent UNSAT | Sum of 32 cold child walls | Charged 32-query portfolio wall* | Max one child wall | Max sampled child-tree RSS |
|:--|--:|--:|--:|--:|--:|
| CryptoMiniSat | 5 / 27 | 9.508 s | 10.326 s | 0.535 s | 102,760,448 B |
| Kissat | 5 / 27 | 6.650 s | 7.468 s | 0.274 s | 130,318,336 B |
| CaDiCaL | 5 / 27 | 10.158 s | 10.977 s | 0.593 s | 136,429,568 B |

\*The portfolio charges the common preflight and 17.6 MB base load once,
each of the 32 derived-file constructions once for that engine, every cold
solver child, and model parsing/point replay. It is a reconstruction from
the frozen sequential three-engine run, not an independent single-engine
rerun. Actual full 96-child Python `main()` wall was 27.316 s, from
14:56:45 to 14:57:12 UTC on the recorded macOS arm64 host; it includes
all 96 children and the in-run setup, but not Python interpreter startup.
The maximum one-child wall is not the portfolio cost. The process-tree RSS
values are 20 ms sampled peaks; Darwin cumulative child `ru_maxrss`
before/after is saved per child but is not a per-child exact peak.

The same fixed T-order first-witness branches were T3,T0,T2,T0 for Q0–Q3,
respectively, for all three engines. Their charged cold Q-prefix wall
(seconds, including common setup) was:

| Engine | Q0 T0..T3 | Q1 T0 | Q2 T0..T2 | Q3 T0 | Q4–Q7 full-negative four-branch wall, respectively |
|:--|--:|--:|--:|--:|:--|
| CryptoMiniSat | 1.511 | 0.381 | 1.088 | 0.320 | 1.573, 1.366, 1.129, 1.111 |
| Kissat | 0.973 | 0.327 | 0.786 | 0.281 | 1.079, 1.026, 1.004, 0.966 |
| CaDiCaL | 1.323 | 0.439 | 1.362 | 0.444 | 1.319, 1.500, 1.321, 1.460 |

For a projected-negative Q all four torsion branches must return UNSAT;
the table charges all four. The positive prefix stops only after a newly
verified SAT witness in fixed T order; all four branches were nevertheless
run and retained. The 5/27 split is **branch** support, while the 4/4
positive and 4/4 negative classifications are **projected Q** decisions.
All per-child CPU, query setup, exact commands and input SHA-256, parse
and certificate wall, output, RSS and both UTC boundaries are in the raw
receipts. No adaptive target label or solver ordering was used.

The 17,630,779-byte base CNF is archived only once in merged #781, SHA-256
`ac6fb610d38023c392308826c88dd5b0c2bbe1347057bd3f29c6a3b9e05d468b`.
Each cold child read a byte-exact derived 17.6 MB file with one positive
unit assumption; its literal, command and derived SHA are saved, and the
archive replay reconstructs that file from the one base. The preregistered
freeze SHA-256 is
`b3b9bbc130e743eb8a803abaa494057909b10cdfc529007d393f1d00c6d7fd82`;
the measured source commit was `cea74c8e0133b972d5cb0e2eff7de2c73feeb308`.
The exact pre-outcome local hash gate and draft PR CI passed before the
first successful smoke or corpus child. A Python 3.9 import-only failure,
a wrong-shim 3.12.13 invocation stopped by the version pin, and a first
3.12.8 tiny-smoke attempt stopped by sandbox-denied Darwin `sysctl` during
psutil child enumeration are all preserved. The final six smokes and
96-query panel used the unchanged frozen runner with exact 3.12.8 outside
the restrictive `sysctl` sandbox; the failed tiny child is not counted as
a corpus outcome. The local host was released before/after the panel for
other timed work.

**Decision:** `ADMITTED_TOY_SOLVER_STAGE` under the frozen cap, with
Kissat having the lowest wall on this one fixed macOS panel. This is an
engineering/solver-admission result, not an index-calculus advance or a
claim of statistical speed superiority. The encoding explicitly lists
rational factor x and O-aware prefix states; it requires 1,195,344
clauses for only n13,m5. The next useful comparison is a separately
verified sparse or algebraic encoding on the same fixed point targets,
then a bounded n19/m6 gate with charged export and independent negative
certification. Dense-CNF size, toy solver wall, and oracle agreement do
not establish n131 tractability, relation rank/yield, end-to-end ECDLP
cost S, or a rho crossover.
