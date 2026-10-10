# N41 C800 combined S3 target path

The combined late root-slot prefetch and projective scalar-check path reduced
the measured one-target online interval by **2.145×** relative to the same
packed S3 index without those changes on eight frozen N41 public points. The
target-level bootstrap 95% interval is **[2.049×, 2.241×]**. All 32 processes
completed rank 800, native group verification, and independent checked-Sage
scalar replay. This is a local algorithm-variant comparison; the host-level
CPU isolation receipt required for a controlled wall-time result is absent.

## Frozen method

The curve is `EC1N41Ckb1hcbf4a37f6e66`, over the polynomial-basis field
defined by `x^41 + x^3 + 1`, with prime subgroup order `549756390943` and
generator `[2056947637384, 1635505394702]`. The low-x factor base contains
65,600 usable subgroup points, folded into 800 signed-Frobenius columns.
Both arms use the same factor-base digest
`cdd82bd0acdb98b53455a97fc5611af7cebdccc87b63d966d9000adb300092a7`,
the same rank 800 setup, the same relation seed `20260928`, 14 preparation
workers, one target worker, and a 32 GiB RSS cap. Each target is supplied as
a public point; fixture construction and reusable setup occur outside its
online timer. The charged online interval includes target query, S3 point
decomposition, relation check, direct descent, and scalar replay.

| Arm | Candidate ID | Target path |
| --- | --- | --- |
| Reference | `IC1N41Ckb1fb65600PDP4rootRCsampleLAgaussTDdirectISO0h161970671e93` | Original packed S3 lookup and affine scalar checks |
| Combined | `IC1N41Ckb1fb65600PDP4rootRCsampleLAgaussTDdirectISO0h700c174acd1e` | Late root-slot prefetch and portable López–Dahab projective scalar checks |

The immutable candidate manifests are under `identities/candidates/`. The
T025 calibration inputs and binary hashes were committed in
`CONTROL_FREEZE.json` before the calibration. Both binaries recovered scalar
`253982070906`, reached rank 800, and returned relation indices
`[32414, 58107, 13391, 43988]`. Their online intervals were 0.244333 ms and
0.089042 ms. The wrapper collecting the reference calibration returned status
1 after writing a valid native report because macOS denied its
`sysctl kern.clockrate` request; the archived stderr and receipt retain that
condition. Calibration was excluded from the eight-point panel.

`CONFIRM_FREEZE.json`, the candidate manifests, and the four runner and audit
sources were committed before the panel. Points T189–T196 follow the parent
freeze's SHA-256 target law. Each point was measured in two fresh processes
per arm, with ABBA or BAAB order. The point is the bootstrap sampling unit.
The local gate required all native and Sage replays, the same first relation
within each four-run block, a geometric ratio of at least 2×, and a target
bootstrap lower endpoint above 1×; the lower endpoint above 2× is recorded
separately.

## Paired result

| Point | Reference online ms, R1/R2 | Combined online ms, R1/R2 | Paired ratio |
| --- | ---: | ---: | ---: |
| T189 | 0.122792 / 0.138291 | 0.063750 / 0.069375 | 1.959× |
| T190 | 0.115875 / 0.114084 | 0.050167 / 0.052333 | 2.244× |
| T191 | 0.179750 / 0.237750 | 0.100375 / 0.088667 | 2.191× |
| T192 | 0.247750 / 0.233541 | 0.118833 / 0.096875 | 2.242× |
| T193 | 0.104625 / 0.123000 | 0.054958 / 0.054917 | 2.065× |
| T194 | 0.184875 / 0.144125 | 0.068708 / 0.069542 | 2.361× |
| T195 | 0.153042 / 0.146625 | 0.067375 / 0.069666 | 2.187× |
| T196 | 0.298042 / 0.297583 | 0.148542 / 0.158209 | 1.943× |

The geometric online ratio is 2.144557×, with target-bootstrap 95% interval
[2.048692×, 2.241089×]. The reference/combined S3-call ratio is 0.850010×;
the combined kernel charges look-ahead calls and does more S3 calls on most
points. The speed gain comes from target-path cost per call and faster scalar
checks. Across the 16 processes per arm, median target-query time was
0.043875 / 0.013937 ms, median charged decomposition time was
0.043188 / 0.030625 ms, and median recovery-check time was
0.035334 / 0.003541 ms (reference / combined). All 32 runs had matching
first relations within their point blocks, and peak RSS was 3,883,171,840
bytes, below the frozen cap.

The 95th percentile absolute A/A log ratio was 0.269 for the reference and
0.176 for the combined arm. Setup time ranged from 0.852–1.128× and
0.873–1.147× its arm median, respectively, failing the frozen 10% host-noise
gate. The host has no exclusive physical CPU partition receipt, so
`controlled_speedup` remains `null` in the audit. The existing independent
N53 C1200 confirmation of these same source variants yielded 1.960× on its
own eight-point panel; curve and base size therefore remain explicit factors
in any comparison.

The complete 141-file archive is under `local/`, including compressed raw
reports, stdout/stderr, per-process resource counters, freeze-linked run IDs,
checked Sage runtime information and scalar replays, the deterministic
20,000-draw bootstrap audit, and `archive_receipt.json`. Reanalysis from the
archived copy produced byte-identical `analysis.json`.

Next, replay this frozen panel on a physical host with audited exclusive CPU
and NUMA allocation, fixed frequency, IRQ routing, and quiet paired runs.
The source-free Modal guest probe in the adjacent experiment did not establish
that host-level isolation, so a Modal timing would retain exploratory status
unless the host provides equivalent evidence.
