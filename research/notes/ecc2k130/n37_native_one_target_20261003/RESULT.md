# n37 descendant-native one-target panel: verified, cost verdict open

Status: **correctness PASS; performance UNDETERMINED**. This is the outcome
of the frozen [protocol](PROTOCOL.md), committed as `97ead2b9` before any
b02 solve. The source curve is ICV1
`icv1-f2m37-tm534059-32aad96b`, with prime subgroup order 230603167;
the IC arm transports it through the certified degree-73 isogeny to the
fixed 42-column descendant-native base. Both arms solved each of 32 new
public Q in its own process, with five paired repetitions per Q and no
shared target-dependent state. The separate source-curve replay checked
all 160 IC/rho pairs, all support and miss decisions, both scalar answers,
and fixture equality only after the producer outputs existed. There were
no failed, timed-out, or incomplete pairs.

| Frozen b02 result | IC fixed K42 + signed 3+3 + 16 shifts | Signed-Frobenius rho, rung 3 |
|:--|--:|--:|
| Verified previously unseen Q | 32/32 | 32/32 |
| Independent paired replays | 160/160 | 160/160 |
| Direct IC hits / shifted recoveries | 27 / 5 | n/a |
| Total IC oracle queries, one run per Q | 39 | n/a |
| Maximum IC queries for one Q | 4 | n/a |
| Median per-Q **recorded explicit point-addition requests**, cold | 1,634,155 | 5,155 |
| Recorded explicit scalar multiplications per process | 197 | 35 |
| Complete-cost `S`, rho ratio, speedup | **unset** | **unset** |

The five shifted Q were indices 0, 6, 10, 25, and 31. The first four used
one shift; Q31 used three. All repeated IC answers, witnesses, complete
misses, base logs, relation rows, and counted work were identical for a
given Q. The relation setup reached rank 42 in 49 target-blind probes;
the signed half table contained 102,391 distinct sums. Per-point counts,
all 160 pair identifiers and file hashes, phases, process observations,
and the five A/A pairs are in [SUMMARY.json](SUMMARY.json).

The addition counters are **not a calibrated common cost unit**. IC counts
batched logical addition requests; rho counts its own explicit additions.
The counts omit scalar-multiplication internals and parts of field,
isogeny, and input work. In particular, rho's public-point subgroup check
performs one scalar multiplication before its online interval that its
existing `charges` counter does not increment; the cold corrected count is
at least 36, rather than the recorded 35. IC's 197 includes its two Q
validation multiplications. The counter gap identifies target-blind
relation search and table construction as a measurement and algorithmic
priority, but dividing the two rows is inadmissible as a speed ratio.

The five exact IC online phases are `target_query`, `target_pdp`,
`target_relation_check`, `target_descent`, and `recovery_check`; their integer
nanoseconds sum to each IC online interval. Rho reports exclusive `walk`,
`collision`, and `recovery_check`, with the same exact-sum property. IC
rebuilds the bridge, base, table, rank-42 relations, and shift points in
each process. Rho uses one Q and a fresh distinguished-point table in
each process. The host envelope is [HOST.json](HOST.json): Apple M4 Pro,
macOS arm64, one Rayon thread, one panel process at a time. This sandbox
could not isolate a CPU core or prove absence of host contention. Five
identical-command A/A pairs had maximum online-time relative spreads of
188.3% for IC and 76.0% for rho. Their timing rows are descriptive only;
no 95% paired runtime claim is eligible. Rho receives Q in an environment
variable extracted by the shell, whereas IC reads and hashes the point
file inside its process, so a complete cold-cost comparison would also
need to charge matched input handling. `S` and speedup remain null.

The hypothesis that the frozen residual policy recovers this 32-Q b02
panel passed. It does not establish a method crossover, a population-wide
success rate, or behavior at n41, n53, m=83, n131, or ECC2K-130. The
research classification is **accounting/correctness**, not an algorithmic
advance. The next controlled decision is to calibrate native field,
batched-addition, scalar-multiplication, and degree-73 transport costs in
one operation unit, and rerun the same one-target protocol on an isolated
host. After that, compare equal-useful-size factor-base policies with
fully charged target-blind relation and table setup. The 49 relation
probes consumed 1,507,265 IC query-addition requests on every Q, so
reducing or reusing this cold setup deserves priority over tuning the
seven b02 shift queries. Reuse is a separately named warm multi-target
question and cannot be credited to the one-target result.

## Reproduction and evidence

The producer/replay source SHA-256 digests are respectively
`b856561d318cdd6ff7083d9c5a96b67901f8844f6c5392d904b0624d0aa57424`,
`f2d0c2eda96864e6d4b39a34830c4c7ed9a7cfb6f18961c8c153d707e012eee3`,
and the instrumented strong-rho source is
`575c1739f37e24d9130cc2fef8ea2c7554c762cf822e52708a13263b2875ccb0`.
The exact `Cargo.lock` bytes are preserved at
`research/notes/ecc2k130/n37_native_m6_mitm_20261002/Cargo.lock`, SHA-256
`b28d3c2d81146a40d00df85f45c6460d9e2bc5c26875307a06cfd15efff99365`.
The input hashes are pinned in [PROTOCOL.md](PROTOCOL.md) and independently
checked by producer/replay. The result archive [RAW.tar.gz](RAW.tar.gz)
contains 1,360 individual raw files (160 paired IC/rho outputs, replay
receipts, process timer logs and stdout logs, plus ten A/A copies), is
1,238,234 bytes, and has SHA-256
`4381af7e409c40cea5b69f102c7d60ee8aa3818f716d8a41b9094c0f9363d78c`.
The 168,349-byte summary has SHA-256
`bef80739aaa96f42f6c34da2ee2b5949bc4dcffdb4906c1dd26f48749d560e10`.
The archive is tracked in Git; a local raw directory is only an extraction
convenience. On a clean checkout, extract it into this note directory and
rerun the native auditor:

```sh
cp research/notes/ecc2k130/n37_native_m6_mitm_20261002/Cargo.lock Cargo.lock
cargo test --release --locked --example koblitz_rho_batch_ks_strong_online
cargo build --release --locked \
  --example n37_native_m6_one_target \
  --example n37_native_m6_one_target_replay \
  --example koblitz_rho_batch_ks_strong_online
sh research/notes/ecc2k130/n37_native_one_target_20261003/run_panel.sh 0 31 0 4
sh research/notes/ecc2k130/n37_native_one_target_20261003/run_aa.sh
cargo run --release --locked --example n37_native_one_target_analyze -- \
  research/notes/ecc2k130/n37_native_one_target_20261003/SUMMARY.json
```

`run_panel.sh` and `run_aa.sh` refuse to overwrite partial outputs. To
verify archived evidence without rerunning solvers, extract the archive
into a clean checkout and run only the final analyzer command:

```sh
tar -xzf research/notes/ecc2k130/n37_native_one_target_20261003/RAW.tar.gz \
  -C research/notes/ecc2k130/n37_native_one_target_20261003
cargo run --release --locked --example n37_native_one_target_analyze -- \
  research/notes/ecc2k130/n37_native_one_target_20261003/SUMMARY.json
```

The analyzer checks every file's
SHA-256 against its independent replay receipt and audits all point,
seed, repetition, answer, counter, and phase invariants. The canonical
[scoreboard](../../../../docs/index-calculus-scoreboard.html) records the
same scoped decision.
