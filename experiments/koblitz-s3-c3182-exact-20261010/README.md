# C3182 exact S3 index at the 8 GiB table boundary

The N53 quotient-root index can use **3182 signed Frobenius orbit columns**
with the existing 8 GiB root table and 29-bit state ordinal. This is the
largest integer column count before both the initial table capacity and the
state ordinal width step upward. The new source preserves the C1200 embedded
target search prefix, so its first 1200-column search follows the frozen
point-seeded state order while querying the larger exact root table.

| Columns | Regular states `53 C²` | State payload | Root slots | Root payload | Combined payload |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 2800 | 415,520,000 | 15.48 GiB | 1,073,741,824 | 8 GiB | 23.48 GiB |
| 3182 | 536,631,572 | 19.99 GiB | 1,073,741,824 | 8 GiB | 27.99 GiB |
| 3183 | 536,968,917 | 20.00 GiB | 2,147,483,648 | 16 GiB | 36.00 GiB |

These are fixed payloads only. The C3182 run has a **32 GiB peak RSS cap**;
the remaining memory must cover the factor base, relation matrix, scatter
buffers, allocator, and runtime. C3182 leaves 239,340 unused 29-bit state
ordinals. The 3183-column case exceeds 2²⁹ states and the next power-of-two
table boundary, so it needs a different representation under the same memory
cap. C3200 falls on that larger-table side as well.

`FEASIBILITY_FREEZE.json` fixes the source and binary SHA-256, C3182 layout,
N53 T025 public point, relation seed, 14 preparation workers, one target
worker, 1200-second limit, and memory envelope before execution.
`run_feasibility.py` runs one target and records complete stdout/stderr,
the raw report, peak child RSS, resource counters, and an auditable receipt.
Its gates require rank 3182 before target input, the exact packed layout,
all five charged online phases, the expected scalar, and group verification.
The target fixture is supplied as a public point. Local wall times from this
host are exploratory until a host-level isolation receipt is available.
