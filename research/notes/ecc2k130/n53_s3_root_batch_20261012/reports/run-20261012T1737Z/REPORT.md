# Valid run record — 2026-10-12 17:37 UTC

The [frozen protocol](../../PROTOCOL.md) requested three new public
hash-to-curve targets, three rank seeds and one matched strong-rho seed per
target, with control, root-batch candidate and rho in rotating order:
27 fresh processes. The frozen source, input and binary identities are in
[`FROZEN.json`](../../FROZEN.json). Execution used an Apple M4 Pro
(`Mac16,7`), macOS 26.6, 48 GiB RAM, Rust 1.93.1 release binaries, with
no core reservation (L0 shared-host timing).

All 27 statuses are `COMPLETED`, with no timeout or observed 16-GiB RSS-cap
breach. Every process kept its command, environment, wall elapsed time,
wrapper maximum RSS, stdout, and stderr at
[`runs/q*/s*/*/`](../../runs/). Native analysis
[`ANALYSIS.json`](../../ANALYSIS.json) verified nine same-target pairs;
all 18 IC arms and nine rho arms verified their scalar. Control and
candidate selected rank traces, base dumps, rank 220/220, 220 relations,
root-table sizes, consumed roots, and target witness agree within each pair.
The exact observations and native medians/intervals are in
[`RESULT_DATA.json`](../../RESULT_DATA.json).

| Measured comparison, same Q (n=9 pairs) | Median ratio | Target-block 95% descriptive interval |
|:--|--:|--:|
| Batch / control, complete cold IC process | 0.737998 | [0.714128, 0.750321] |
| Batch / control, target-online IC phase | 0.611643 | [0.568613, 0.903106] |
| Control / strong rho, complete cold process | 78.145833 | [44.219726, 92.601592] |
| Batch / strong rho, complete cold process | 52.585601 | [33.179008, 64.667877] |
| Control / strong rho, target online | 0.489747 | [0.167802, 0.763394] |
| Batch / strong rho, target online | 0.377771 | [0.102635, 0.442294] |

The chart [paired ratios](paired_ratios.svg) plots all nine recorded
candidate/control ratios and native medians; [rho reference](rho_reference.svg)
plots the same-point rho ratios. Both axes are unitless ratios of
in-process milliseconds. The displayed intervals are native deterministic
target-block bootstrap summaries, not host-isolated performance confidence
claims. [This PDF](RUN_REPORT.pdf) embeds both figures. Source-level
rank-root sums across the nine cells are 80,790,012 versus 763,167 inversion
primitives, 1,211,850,180 versus 902,063,394 field multiplication calls,
and 242,370,036 versus 39,684,684 field squaring calls (control versus
batch). The method-wide common operation unit `S` is unmeasured.

This run report records observations. Gate interpretation and unresolved
common-unit/host checks are in [`RESULT.md`](../../RESULT.md). The preceding
sandboxed wrapper attempt is a separate infrastructure failure under
[`runs_sandbox_wrapper_failure_20261012T1735Z/`](../../runs_sandbox_wrapper_failure_20261012T1735Z/).
