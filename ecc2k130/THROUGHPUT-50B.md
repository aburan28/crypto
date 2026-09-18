# 50 B/s and fruitless-cycle boundary

Target: **50 billion iterations/s** on one RTX PRO 6000, with no unpriced
fruitless-cycle behavior. Result: **not reached**.

## Boundary

Unit for primitive rows is complete operations/s. Unit for the mixed row is
effective rho iterations/s after the published r+h randomness factor.

The one-product lambda-halving kernel now issues 14 CLMAD instructions per
half. At the measured 1.684 lanes/SM-clock, 188 SMs and 2.4 GHz, its absolute
CLMAD-only ceiling is 54.27 B/s with every ALU instruction, load and dependency
free. A 50 B/s target would require 92.1% of that absolute issue ceiling while
also executing the half-trace, square root, trace test and state traffic.

## Single table

| variant | measured/model B/s | / 50 | correct | cycle status | class |
|---|---:|---:|---|---|---|
| selected affine table-add rho step | 20.134 | 0.403 | 300/300 | exact add-tag 2/4-cycle escape | reference |
| one-product lambda halving, 768 threads | 28.953 | 0.579 | 512/512 | permutation, not rho | engineering |
| one-product lambda halving, 1024 threads | **29.053** | **0.581** | **512/512** | permutation, not rho | **engineering; target not met** |
| explicit dual-chain ILP | 28.823 | 0.576 | 512/512 | permutation, not rho | engineering, rejected |
| binary-tensor square-root map | 1.746 | 0.035 | 512/512 | permutation, not rho | engineering, rejected |
| CLMAD-only halving ceiling | 54.273 | 1.085 | extrapolation | excludes map semantics | hardware bound |
| ideal r=128 mixed optimum | 22.257 raw / **21.112 effective** | 0.422 effective | model | cycle guard unpriced | upper model |
| 50 B/s target | 50.000 | 1.000 | required | zero undetected fruitless cycles required | boundary |

The mixed model optimizes the halving probability at 0.31071. Its operation
harmonic mean is 22.257 B/s and its r+h randomness factor is 1.05424, leaving
21.112 effective B/s before lambda/affine conversion, GPU branch dispatch,
cycle detection or restart. Thus a 50 B/s primitive would not imply a 50 B/s
rho walk.

## Fruitless cycles

Literal “no cycles” is impossible for Pollard rho: every deterministic map on
a finite set is eventually periodic, and collision finding relies on that
fact. The enforceable requirement is **no undetected fruitless short cycle and
no trail allowed to spin indefinitely without a distinguished point**.

The existing additive table walk rejects exact inverse-tag 2-cycles and
four-cycles using the last three tags. A mixed halving walk has not been
enabled. Before it can count as a result it must add, price and test:

1. exact coefficient replay for both `Y+M_j` and `[2^-1]Y`;
2. the existing addition 2/4-cycle escape;
3. exact state-cycle detection with restart for cycles involving halving;
4. a maximum-trail guard as a backstop;
5. exhaustive small-curve functional-graph tests showing every detected
   non-DP cycle restarts and zero such cycles are counted as relations;
6. measured cycle-detection and representation-conversion costs in the same
   throughput row.

Until those gates pass, the 29.053 B/s number is a cheaper primitive only.
There is no 50 B/s rho claim.

Binary tensor cores were also measured, not assumed. A direct 131×131
polynomial square-root matrix using exact XOR-popcount BMMA passes all 512
points but runs at 1.746 B/s: 68 warp-level matrix operations plus fragment
shared-memory traffic per 32 points cost far more than the scalar permutation
network. It is retained and rejected.

Frozen evidence:
[`benchmarks/throughput-50b-halving`](benchmarks/throughput-50b-halving/).
