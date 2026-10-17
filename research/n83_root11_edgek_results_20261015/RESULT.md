# N83 root11 capacity-edge factor-base screen

K3450/root11/wide3 gave the lowest mirrored target-online mean, **20.426 ms**
on 1,000 fresh public subgroup points. Increasing the factor base to K3500
or K3550 increased the measured mean despite adding 8,300 or 16,600 usable
points. K3550 was 2.171 ms slower than K3450 with a contiguous-block 99%
interval of [+0.182, +4.252] ms. The K3500 minus K3450 interval was
[-0.512, +3.164] ms, so this panel records `BOUNDED_UNKNOWN` for a single
winner. Both larger K means are above K3450, so the prospectively frozen
`../n83_edge_confirmation_20261015/PROTOCOL.md` admits no edge challenger.

| Root11/wide3 pipeline | Usable base points | Folded columns | Exact root keys | Mirrored target-online mean | Mean minus K3450, 99% interval | Peak RSS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| K3450 | 572,700 | 3,450 | 987,887,478 | 20.426 ms | reference | 11.504 GiB |
| K3500 | 581,000 | 3,500 | 1,016,729,255 | 21.827 ms | +1.401 ms [-0.512, +3.164] | 11.656 GiB |
| K3550 | 589,300 | 3,550 | 1,045,986,026 | 22.597 ms | +2.171 ms [+0.182, +4.252] | 11.813 GiB |

The six complete pipeline runs used list/reverse-list order, a 14 GiB RSS
guard, four rank threads and one target thread on one cloud allocation.
Each run built a fresh base and exact root index, reached full rank, solved
its 1,000 public points, and passed independent base, four-point relation and
scalar replay. All runs recorded zero new CPU-throttled periods and OOM kills.
Setup took 470–536 seconds per run; full process wall time was 522–591
seconds. The table reports target-online intervals after reusable setup;
`measurements.json.gz` retains exclusive target phases and all six setup,
rank, final LA, memory, probe and status records. This is a secondary
factor-base selection screen under the stated cloud-instance envelope.

The public target SHA-256 is
`c75cbb92ef1a23057272cf3481cc649b12495135c75affd1c2610205a42b7aef`;
the answer sidecar SHA-256 is
`ca95d57ff091b8fc7de8aa2fd22120cadac94a7c7a697e3c757aebc53115a159`.
The producer source SHA-256 is
`f2745227d95e36ef45c73c5104a23ada6cebb0794405368581d2be17044845dc`
at commit `10c7c80d9489ef4fd1f302c9063c7452444aeaee`.
[Cloud run 38794757719](https://github.com/aburan28/crypto/actions/runs/38794757719)
uploaded artifact 12005586706. The archived ZIP SHA-256 is
`a7bbcefaff56322e234db96ea799771941dab020252e0adc89a6c66190ba62e3`.
Cloud and separate local replay produced the same analysis SHA-256,
`f50192442355a64fd1c4cfd07147f5b35f8aef4ae77d783f1c7614035992c3eb`.
The archived ZIP and expanded raw receipts live on the second local volume
under `/Volumes/SSD990-2/n83-archive-20261015/raw/n83-root11-edgek-38794757719`.
