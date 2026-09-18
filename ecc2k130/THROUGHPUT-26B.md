# Pursuing 26 B/s on one RTX PRO 6000

Target: **26 billion complete scalar updates/s**. Result: **not reached**.
The verified reference remains 20.134 B/s.

## Boundary

Unit: billions of complete scalar updates per second. The current table walk
performs one affine addition per update. Its derived practical boundary on
this GPU is 22--25 B/s, with the upper endpoint already assuming nearly ideal
pipe utilization. Therefore 26 B/s is at least 1.04 times the upper boundary
and cannot be reached by selector, scheduling or cache tuning of the same
iteration.

Success required a verified median above 26, 300/300 replayed reports, zero
drops and complete accounting. An alternative iteration could move the
boundary, but its effective collision work must also be priced.

## Single table

| variant | measured B/s | / 26 | / 20.134 reference | correct | class |
|---|---:|---:|---:|---|---|
| selected B17 table walk | **20.134** | 0.774 | 1.000 | 300/300, 0 dropped | reference |
| fused 16-bit weight/phase LUT | 20.055 | 0.771 | 0.996 | 300/300; 4096/4096 selector points | engineering, rejected |
| exact DP4A phase selector | 19.515 | 0.751 | 0.969 | 300/300; 4096/4096 selector points | engineering, rejected |
| two-product point halving | 16.053 | 0.617 | 0.797 | 512/512 subgroup halves | alternative primitive; not a rho map alone |
| one-addition practical boundary | 22--25 | 0.846--0.962 | 1.09--1.24 | derived | boundary |
| **26 B/s target** | **26.000** | **1.000** | **1.291** | required | above boundary |

## DP4A result

`TABLE_DP4A_PHASE=1` computes the exact Frobenius phase with 33 unsigned
four-byte dot products. Each four-bit input nibble expands to four 0/1 bytes,
which are dotted with packed logarithm weights. It replaces 17 random shared
byte lookups and shrinks the phase table from 4,352 to 132 bytes.

The selector is exact, but the arithmetic costs more than the conflicted LUT:
three samples were 19.510, 19.530 and 19.515 B/s. The 19.515 median is a 3.1%
regression. The device used 102 registers versus 98 for the selected LUT
kernel. It remains an opt-in negative result.

The second selector candidate stores `log_sum + 4096*popcount` in each 16-bit
entry, deriving phase and Hamming weight from the same 17 loads and deleting
five quarter-rate POPCs. It retains 98 registers and is exact, but the wider
shared loads enlarge the table to 54,132 bytes and measure 20.055 B/s, 0.4%
below the selected byte-LUT row.

Point halving did change the product count, but its linear solves and subgroup
selection measured only 16.053 B/s, and halving alone is a permutation rather
than a rho iteration. Thus neither measured alternative moves effective
throughput past the one-addition boundary.

Frozen build, correctness and benchmark logs are in
[`benchmarks/throughput-26b-attempt`](benchmarks/throughput-26b-attempt/).
