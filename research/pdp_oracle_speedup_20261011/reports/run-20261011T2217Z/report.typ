#set page(margin: 1.8cm, numbering: "1")
#set text(size: 9.5pt)
#set heading(numbering: "1.")
#show raw: set text(size: 8pt)

#align(center)[
  #text(size: 15pt, weight: "bold")[PS-1 run report: algebra spec, bit-degree and root-test formulations] \
  #text(size: 10pt)[PDP oracle speedup program, `research/pdp_oracle_speedup_20261011/ps1_algebra/` --- run of 2026-10-11, 22:17 UTC]
]

= What ran

`ps1_algebra/run_ps1.sh` with `examples/pdp_speedup` (subcommands `degrees`, `torsion`, `agree`, `bs-check`, `sample`, `corpus-check`, `gates`), Conductor task T-4, on the PS-0 host (Intel Xeon 2.80 GHz KVM guest, AVX-512, no GPU). No timing is taken; every figure is a count or an exact comparison. Raw records: `ps1_algebra/raw/*.jsonl`.

= Observations

*Bit-degree (measured, $l = 3..8$, six targets each).* Joint algebraic degree of the output bits of $c_0..c_4$ in the bits of $(x_1, x_2)$: 2, 4, 2, 4, 2 at every $l$; in the bits of $x_2$ alone with $x_1$ fixed: 1, 2, 1, 2, 1.

*Exact agreement (summation-polynomial roots against the group-law reference).*

#table(
  columns: (auto, auto, auto, auto),
  align: (left, right, right, right),
  table.header[*check*][*cases*][*items*][*disagreements*],
  [exhaustive, fields of degree 5--9, all $l <= n$, all rational targets], [35], [4168 targets, 15,628,190 roots], [0 at non-degenerate targets; 1300 explained extra roots at degenerate targets; 0 unexplained],
  [bit-sliced 64-lane pipeline against exact reference, same range], [35], [2834 targets, 1,259,432 pairs], [0 (10,814 zero-denominator pairs routed to the exact path, not compared)],
  [sampled, $(n, l) in {(11,3),(11,4),(11,5),(13,4),(13,5),(13,6)}$], [6], [$6 times 10^6$ random pairs], [0],
  [PS-0 corpus root sets, $l = 8..12$], [5], [400 targets], [0],
  [translation by the 2-torsion point is $x |-> 1/x$, degrees 5--12 exhaustive, 24--36 at 200,000 samples], [13], [about 1.0 million points], [0],
)

#figure(
  image("gates_per_task.svg", width: 90%),
  caption: [Gates per evaluation of the root test, by formulation and rung (source `ps1_algebra/raw/gates.jsonl`; fused ternary-logic count). The baseline formulation is derived from PS-0's counted operations; the group-law task is counted by running the same code on a gate-counting lane. Log scale.],
)

#table(
  columns: (auto, auto, auto, auto, auto, auto),
  align: (right, right, right, right, right, right),
  table.header[*l*][*n*][*group-law task, fused*][*unfused*][*baseline (derived)*][*ratio*],
  [8], [24], [7,410], [12,766], [70,695], [9.5],
  [9], [27], [9,137], [15,961], [94,729], [10.4],
  [10], [30], [10,354], [18,869], [118,099], [11.4],
  [11], [33], [12,349], [22,653], [151,371], [12.3],
  [12], [36], [14,608], [27,011], [190,275], [13.0],
)

= Missing or not done

No `identity.certificate/v1` record exists: the schema is fixed to $F_p$ with $p = 2^(61) - 1$ and cannot express a characteristic-2 identity. The baseline formulation's count is derived from counted operations, not a built program. The trace-based and resolvent-cubic candidates were not attempted. The sampled comparison at $n = 11, 13$ is against the generic summation-polynomial enumerator, not against `decompose`. Wall time was not measured; the group-law test is a different algorithm, valid for targets outside $x_R in V union V^(-1) union {0}$.
