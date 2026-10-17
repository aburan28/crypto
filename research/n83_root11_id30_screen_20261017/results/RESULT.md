# N83 five-byte ID30 high-K screen

The five-byte ID30 root index constructed the K4096/root11 factor base with
679,936 usable points, 4,096 folded columns, full rank, and independent replay
of all 128 four-point target relations and scalar answers. Its target-online
geometric mean was 22.466029 ms. The frozen K3072 control measured
16.405336 ms on the same 128 public points, giving a paired geometric
K3072/K4096 ratio of 0.730229 (descriptive 95% interval
[0.544154, 0.979931]); K4096 was faster on 56 of 128 points. The expanded
base reduced mean lookup probes to 36,347.24 from 58,734.66 at K3072, but
its larger retained state and root table outweighed that reduction in this
local target-online measurement.

| Frozen arm | Usable points | Folded columns | Producer status | Target-online geometric mean | Cold process | State plus root index | Child maximum RSS receipt |
| --- | ---: | ---: | --- | ---: | ---: | ---: | ---: |
| K3072 control | 509,952 | 3,072 | independently replayed prior row | 16.405336 ms | 491.832 s | see control artifact | see control artifact |
| K4096 five-byte ID30 | 679,936 | 4,096 | `EXIT_0`; 128/128 replay | 22.466029 ms | 1,500.010 s | 17.133 GiB | 14,980,890,624 bytes |
| K4608 five-byte ID30 | unknown | unknown | `TIMEOUT` at 1,800.734 s during construction | unknown | capped | no completed allocation receipt | cumulative child maximum 20,909,326,336 bytes |
| K5086 five-byte ID30 | unknown | unknown | `NOT_RUN_AFTER_STOP_RULE` | unknown | unknown | unknown | unknown |

The K4608 producer emitted no completed base dump, summary, or target rows
before its predeclared cap. Its status is a censored setup measurement; the
20.9 GB `RUSAGE_CHILDREN` value is cumulative across the K4096 and K4608
children and is an upper-bound observation for K4608. The stop rule then
left K5086 unmeasured. The K4096 row's 10 GiB root table and approximately
7.133 GiB retained states explain the 17.133 GiB allocation figure; the
maximum resident memory receipt is distinct from that allocation count.

This was an ordinary local diagnostic with two sequential fresh processes,
not a hardware-qualified wall-time comparison. The same 128-point input and
answer sidecar were frozen by plan SHA-256
`ce79d7c6ebadb7a5d38a4d4b9609dabda3e3a224c4a7409862533159f0f07938`.
The Rust source SHA-256 was
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`;
the release binary SHA-256 was
`faaac22aa6228c8d486875d6947fa654417d1122191fbe0eb9bde59bb45eb98f`.
The exact [summary](summary.json), [K4096 replay](k4096_analysis.json),
[K4608 capped receipt](k4608_timeout_receipt.json), and
[plan receipt](plan_receipt.json) are versioned here. Full raw commands,
stdout, target rows, base dump and stderr remain under
`/Volumes/SSD990/crypto/artifacts/n83_root11_id30_screen_20261017`.

The next factor-base screen uses the frozen K3088, K3104, K3114 and K3136
grid around the K3072 root-table capacity edge. Separately, the checked
K3072 target-cancellation combination tests whether two cooperative workers
lower target-online time on the current best measured local base. The
four-byte ID30 K4096 layout tests whether reducing root-table allocation
changes the larger-base result.
