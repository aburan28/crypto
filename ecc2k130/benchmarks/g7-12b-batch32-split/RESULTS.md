# Larger split-batch result on G7

Neither larger logical batch advances. The selected B16 runtime remains
unchanged, and the 12 B/s goal remains open.

| Variant | Logical workers | Physical slots/thread | Median B updates/s | Paired ratio | 95% CI | Rate / 12 B/s | Decision |
|---|---:|---:|---:|---:|---:|---:|---|
| selected B16 | 524,286 | 8 | 6.274186 | 1.000000 | [1.000000, 1.000000] | 0.522849 | reference |
| split B24 | 349,524 | 12 | 6.121832 | 0.971864 | [0.961912, 0.981919] | 0.510153 | regression |
| split B32 | 262,143 | 16 | 5.994967 | 0.948253 | [0.930625, 0.966215] | 0.499581 | regression |

Every one of the nine retained samples advances 8,388,576 scalar walks by
4,096 steps, or 34,359,607,296 complete updates, on the same AWS g7.2xlarge /
RTX PRO 4500 at 165 W. Process observations found no foreign GPU process, and
all runs report the exact update count and zero drops. Paired intervals use
Student t on three log ratios.

The B16 isolated build matches every selected native function. B24 and B32
retain 80 registers, an 8-byte spill frame, 29,568 shared bytes per block and
three resident blocks per SM. B32 also retains the selected walk's 5,952
static instruction sites and opcode counts. The larger batches nevertheless
increase median memory activity from about 12% selected to 36% for B24 and
50% for B32; sustained SM clock falls from about 1,942 MHz to 1,852 and
1,755 MHz at the same power cap. Caching only four X slots covers one third of
a B24 physical batch and one quarter of B32, so extra coordinate traffic costs
more than the saved root-inversion work.

All built-in suites pass. Exact cross-batch normalization covers ragged and
block-aligned populations, 86,600 matching DP records and 96 CPU replays.
Both candidates resume bit-for-bit within their own checkpoint format and
reject B16 checkpoints. Six storage/inverse/actual-walk units and six complete
walk sanitizer runs pass with zero errors and zero drops.

Generic work remains `sqrt(n/262)`, ratio 1. Full-DLP S is null.

[Frozen comparison](comparison.json) · [Validation audit](validation-audit.json) · [Native audit](native-audit.json) · [Protocol](README.md)
