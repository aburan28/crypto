#set page(paper: "a4", margin: 16mm)
#set text(size: 10pt)
= Independent Linux isogeny fiber replay
*GitHub Actions run 38831567671 · 2026-10-16 00:50–00:54 UTC.*
Commit `de340fa2161088a5188b841ec82278f5067f8aed` built the frozen
native runner and checker on an Ubuntu 24.04 Azure hosted VM with four
AMD EPYC 9V74 vCPUs. Every one of 252 cold children exited zero. The
checker replayed 5,355 scalar outputs and 37,061 relation witnesses;
the host analyzer matched all 252 counted arithmetic vectors to the
committed Darwin run. The raw 1,330,324-byte archive and its verified
1,533-file manifest are in `linux_host_v1/`. Full source links and
SHA-256 hashes are in `REPORT.md`.

At 64 public targets, median whole-child process CPU ratios over three
paired rounds were full-image/rho 68.750, 56.331, 59.050 at subgroup
orders 421, 433, 457. Full-image/BSGS was 135.798, 109.094, 114.506;
rho/rho-A/A was 0.994, 1.006, 0.994. Whiskers show each three-round
minimum and maximum. The hosted VM has no reserved physical core, so
this is descriptive cross-host CPU evidence, not admitted isolated wall
performance or complete portable addition-equivalent cost.

#image("GRAPH.svg", width: 100%)
