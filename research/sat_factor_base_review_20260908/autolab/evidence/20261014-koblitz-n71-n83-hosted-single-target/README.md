# n=71 and n=83 compact orbit: hosted isolated single-target reruns (2026-10-14)

Status **`PENDING_INDEPENDENT_VALIDATION`**. These are independent-host reruns
of the two single-target compact-orbit autolab beats on hosted GitHub runners
(AMD EPYC 7763, 4 vCPUs, Linux, `ubuntu-24.04`), dispatched with
`workflow_dispatch` from `fix/single-target-panel-point-encoding`
(run [38585913645](https://github.com/aburan28/crypto/actions/runs/38585913645)
for n=71 and
[38597929540](https://github.com/aburan28/crypto/actions/runs/38597929540)
for n=83) under `tools/isolated_bench.py reserve` with every producer pinned
to the reserved CPU. Each run replays its own records on the same runner.
These are independent-host reruns of the measurements, not
independent-author replays, and they promote no ledger row.

The beats are one-target workloads on public hash-to-curve points with no
known scalar: verification is by cross-arm agreement (IC and rho recover the
same scalar on the same target), plus out-of-process Python certificate
replay. The comparison arm is strong rho rung 3
(`koblitz_rho_batch_ks_strong_online`), 32 lanes, empty DP table per
process.

## Why these runs needed the panel fix

The local macOS runs of the same beats finished `VERIFICATION_FAILURE`
despite both arms recovering identical scalars in 16/16 paired workloads:
`rho_ok` compared rho's `published_q` (decimal strings from the wide u128
producer) against the panel's Python-int targets with raw equality, and the
IC relation replay fed string base-dump coordinates into integer field
arithmetic. PR #1827 normalizes point encodings; both hosted runs executed
with that fix and reach `PENDING_INDEPENDENT_VALIDATION`.

## Headline (hosted, isolated)

| Beat | Run id | K | W | online speedup rho/IC (median, 95% CI bootstrap) | cold IC/rho | IC online median ms | rho online median ms |
|---|---|---:|---:|---|---:|---:|---:|
| n=71 a=0 | `gha-38585913645-1` | 600 (fixed) | 16 | **7.29x** (4.64–13.38) | 154.4 | 518.6 | 3764.5 |
| n=83 a=1 | `gha-38597929540-1` | 600 (fixed) | 16 | **12.29x** (3.26–32.21) | 133.7 | 557.5 | 5442.8 |

`verify` PASS for both runs (256 files, 0 mismatches). `isolation.jsonl`:
1 reservation sample each, **0 contended**. Claim-checks: 16/16 PASS per
run. K is fixed at 600 (the pick of the local macOS tune
`20261014T032514Z-ad1a92c1fb` for n=71 and `20261014T122948Z-289481dc5c`
for n=83), so `k_source: fixed` and neither hosted run has a tune table.

## Local (contended) comparison runs

| Beat | Run id | Host | online speedup rho/IC (median, 95% CI) | notes |
|---|---|---|---|---|
| n=71 a=0 | `20261014T032514Z-ad1a92c1fb` | Apple M4 Pro, macOS, unpinned shared host | 5.6x (3.6–10.3) | 16/16 claims, all replays pass; contended host, wall secondary |
| n=83 a=1 | `20261014T122948Z-289481dc5c` | Apple M4 Pro, macOS, unpinned shared host | 9.46x (2.53–34.44) | 16/16 claims, all replays pass; host load 15–127 during run |

Hosted and local medians agree within their bootstrap intervals; the hosted
numbers are the ones earned under isolation.

## Scope and non-claims

- Public synthetic hash-to-curve targets only; no key recovery, no
  asymptotic sub-rho claim, no deployed-curve-security statement.
- The IC's reusable setup (about 500 s in process for n=71, mostly rank
  PDP) is outside the online interval and reported beside it: the cold
  IC/rho ratio is ~154x (n=71) and ~134x (n=83) with setup paid once per
  process. Amortizing the setup across targets is a separate question.
- Wall time on one host class (AMD EPYC 7763, Linux) plus one contended
  Apple M4 Pro; no claim for other hardware classes.
- Out-of-process Python replay on the same runner; independent-host
  *author* replay and a Bernstein–Lange reference remain outstanding.
- Both arms ~sqrt(L*r/n) at optimal K; the online gap is a constant-factor
  result at these sizes, not an asymptotic crossover.

## Files

- `autolab_runs/gha-38585913645-1/` — n=71 hosted run (logs, artifacts,
  claims, receipts, state).
- `autolab_runs/gha-38597929540-1/` — n=83 hosted run.
- `verify.json`, `verify-n83.json` — `boundary_autolab.py verify` receipts
  (PASS, 256 files each).
- `isolation.jsonl`, `isolation-n83.jsonl` — `tools/isolated_bench.py`
  reservation records (0 contended).
- `lscpu.txt`, `lscpu-n83.txt` — host manifests.
- `figures/` — recorded per-workload online walls
  ([SVG](figures/n71-n83-hosted-online-ms.svg),
  [PNG 150 dpi](figures/n71-n83-hosted-online-ms.png),
  [data](figures/data.json), [source](figures/plot.py)).
- `REPORT.typ`, `REPORT.pdf` — the report rendering with the figure.

## Reproduce

```sh
gh workflow run compact-orbit-n61-single-target-isolated.yml \
  -f beat=koblitz.compact_orbit.n71_single_target \
  -f workloads=16 -f tune_workloads=4 -f k=600 -f timeout_minutes=350
```

(same with `beat=koblitz.compact_orbit.n83_single_target` for n=83). The
workflow builds the producers from `Cargo.lock`, reserves a core with
`tools/isolated_bench.py`, and runs
`boundary_autolab.py launch-single` end to end, then `verify`.
