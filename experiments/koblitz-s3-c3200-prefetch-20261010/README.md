# C3200 mixed quotient-root target prefetch

The C3200 target extractor now computes up to eight states' ordered S3 roots,
prefetches the corresponding mixed quotient-table home slots, then performs
the same exact lookup and group-lift sequence. The first verified relation
and all scalar checks remain in the target online interval. This changes the
target PDP memory-access schedule while retaining the C3200 factor base,
root-table representation, preparation rank, and C1200-first traversal.

The [frozen T025 control](CONTROL_FREEZE.json) binds source SHA-256
`56ed0e2167c9b5b1cc389b4d95ed7887f65bf17273d68911b9924ee4167008ac`,
release binary SHA-256
`54644a12e05ab3212dbad3dbf723f66dc913d1fb66bdaec135128570abb3dfae`,
the existing C3200 projective report, public point, runner, 14 preparation
workers, one target worker, 1,200-second cap, and 32 GiB observed-RSS cap.
Its pass criteria include the same first relation witness, 339,200 usable
base points, rank 3,200, exact mixed quotient table, the five exclusive
online phases, and recovered scalar. All 17 release example tests passed.

The control passed: rank 3,200, the same first relation witness
`[55035, 67241, 57152, 189185]`, independent Sage scalar replay, and
28,387,180,544-byte peak RSS. Its T025 target online interval was 1.450375 ms,
versus 3.036625 ms for the earlier projective control; target PDP was
1.353333 ms versus 2.939542 ms. The prefetch batch completed 512 S3 calls
versus 498 in the reference. These single-run timings are diagnostics on
an unisolated host.

The [frozen paired screen](SCREEN_FREEZE.json) compares this source with
the matched-rebuild C3200 projective reference on the eight previously used
T173–T180 public points, in alternating arm order. Both arms use the same
14 preparation workers, one target worker, 3,200-column rank, and 32 GiB
RSS cap. The primary screen statistic is the geometric mean of paired
single-target online ratios; every run must verify and fit the memory cap.
The screen advances only at 2×. Local CPU ratios remain exploratory until
host-level isolation and noise checks pass.

## Paired screen result

All 16 runs verified the supplied public points, stayed below the 32 GiB
observed-RSS cap (maximum 30,934,007,808 bytes), and passed independent Sage
scalar replay. Both arms found the same first relation witness on every
point. The paired online geometric mean was **1.167×** in favor of prefetch,
with target-resampling bootstrap 95% interval **[0.401×, 3.386×]**; the
predeclared 2× development gate did not pass. The geometric mean of baseline
over candidate S3-call counts was 0.946× because the candidate completes
its final eight-state batch before searching the prepared lookups.

| Point | Reference online ms | Prefetch online ms | Reference / prefetch | S3 calls reference / prefetch |
| --- | ---: | ---: | ---: | ---: |
| T173 | 0.097 | 0.680 | 0.14× | 14 / 16 |
| T174 | 1.307 | 0.292 | 4.47× | 80 / 80 |
| T175 | 0.311 | 0.920 | 0.34× | 142 / 144 |
| T176 | 0.556 | 0.090 | 6.17× | 42 / 48 |
| T177 | 0.204 | 0.775 | 0.26× | 242 / 256 |
| T178 | 1.254 | 0.407 | 3.08× | 108 / 112 |
| T179 | 0.648 | 1.850 | 0.35× | 408 / 416 |
| T180 | 1.016 | 0.112 | 9.10× | 122 / 128 |

The [post hoc order audit](local/pair/posthoc_order_audit.json) found a
systematic host effect: the second run was slower on all eight points, with
4–5 major page faults over the process, while every first run recorded zero.
Second-over-first online time had a 4.518× geometric mean. The process-level
fault counter does not localize those faults to the target interval, but
target validation and target PDP times show matching spikes. Reusable setup
times also failed the ±10% diagnostic in both arms. The alternating schedule
balances arm order, yet these observations prevent a controlled CPU speedup
claim from this host.

The [analysis](local/pair/analysis.json), [completion](local/pair/completion.json),
[Sage replay](local/pair/independent_sage_replay.json),
[Sage runtime receipt](local/pair/sage_runtime_info.json), and
[complete raw archive](local/pair/raw.tar.gz) retain the results and failures.
`archive_pair.py` verifies every extracted file against the source tree when
creating the archive. The next performance experiment should use an audited
isolated host and focus on the stronger C1800 combined candidate; this
C3200 prefetch variant remains a correctness-validated diagnostic.
