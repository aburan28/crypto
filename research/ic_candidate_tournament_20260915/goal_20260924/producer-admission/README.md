# Optimized producer admission evidence

The three restored optimized producers pass correctness, canonical identity and
exclusive instruction-accounting checks. This is an admission result, not a
comparative performance result or an improvement round. No incumbent has been
selected and no speedup or globally best implementation is claimed.

Implementation: [PR 733](https://github.com/aburan28/crypto/pull/733).
Tested implementation head: `6b5a402558c2c883baa52c2ef585f99886108167`.
[Final Linux run](https://github.com/aburan28/crypto/actions/runs/36106499204)
uses Rust 1.94.1, static musl and Valgrind 3.22.0. Each process has one pinned
CPU, an 8 GiB address-space cap and a 180-second limit. Source jobs run on
separate CI hosts, so their timings must not be used as a paired comparison.
The [protocol](../../producer/PROTOCOL.md) was committed before these controls.

| Restored source | Cells | Verified native/profile pairs | Derived source-manifest SHA-256 |
|---|---|---:|---|
| `both` | 13a0, 17a1, 23a1 | 9 / 9 | `fdb0704707b7b5125a30123ab45d8d0a90d0c6da2df90996e17bc2029eab0f53` |
| `scaled` | 13a0, 23a1, 37a0, 43a1, 61a1 | 15 / 15 | `3fa7511e1cb1a83ad23d869e4a6bcc4fad3e6fcc8c4d4fd8cb89f98674c85ce3` |
| `pairinv` | 13a0, 23a1, 37a0, 43a1, 61a1 | 15 / 15 | `205e4b2d58e883b3b84b8b36634d8a71f334d770bfa0aa619c5c78e8042fe3cd` |

There is one fixed public target per cell and three fresh process repetitions.
Repeated processes are not fresh targets, independent yield samples or held-out
confirmation. The 39 pairs cover 13 source/cell combinations and six distinct
curve cells. They do not satisfy the goal's 60 fresh confirmation-target gate.

The transported evidence was audited independently of the running solver using
[audit_transport.py](audit_transport.py). It checks the complete source inventory,
measured executable, frozen evaluator, protocol, every receipt's file hashes,
base census, curve/candidate/workload identities, query chronology, rank events,
matrix digest, group relations, factor logs and recovered target scalars. It
reconstructs every canonical run from the raw native report and Callgrind files;
all 39 records agree and all exclusive phase sums close against the profiler's
whole-process instruction checksum. The reconstruction receipts are in
[audits](audits/).

Additional validation on the tested head: 83 Python harness tests pass on Linux;
optimized release controls pass for all three sources (9 tiny-IC and 4 worker
tests for `both`, 10 tiny-IC and 4 worker tests for each later source). The
[existing autolab CI run](https://github.com/aburan28/crypto/actions/runs/36106499266)
also passed. Its 34 reports, 32 base inventories and four native/profile pairs
were independently rechecked after download. Rust lint and workflow syntax pass.

## Retained failures and corrections

[history.json](history.json) identifies every producer CI attempt. The raw archive
retains the first two runs, their failures and logs alongside the passing run.

- Run `36105578831` exposed a stale-row iteration bug in the new Python rank
  trajectory audit. Replacing a reduced row did not replace the list iterator;
  a later pivot used an old coefficient. A real degree-23 vector now guards this
  case. The existing independent final-rank checker correctly accepted the
  complete scalar-field matrix.
- The archived `both` release test still multiplied the descent scalar by the
  cofactor although its columns already represent subgroup points. Its group
  witness and recovered scalar checks passed; the obsolete scalar assertion
  failed. The later `scaled` source already corrected this convention. The
  derivative applies the same test-only correction and retains all group and
  scalar checks. Original archived sources remain unchanged.
- Run `36105859677` passed all `scaled` and `pairinv` checks but retained the
  `both` assertion failure. It is not a passing overall qualification run.
- The initial artifact uploads omitted `.cargo/config.toml` because hidden files
  were excluded. Those incomplete transports remain unchanged and are explicitly
  identified in the history. The final run includes the complete source tree,
  measured executable and frozen evaluator.

The source review also corrected the inherited winner-summary source hash,
required the actual `tiny_gauss` dispatch, removed biased nonzero-scalar draws,
required nonempty IC descent witnesses and added compiled source tags after a
locally reproduced Cargo-cache mix-up. These are accounting/admission corrections,
not measured algorithmic gains.

## Restore and replay

The durable archive is
[`ic-producer-admission-20260925.tar.zst`](../../evidence/ic-producer-admission-20260925.tar.zst).
Its SHA-256, size and file count are in the existing
[evidence manifest](../../evidence/manifest.json). The archive includes raw
profiles, stdout/stderr, source snapshots, executable, frozen evaluator,
canonical records, logs, failure history and audit receipts. Build caches and
Python bytecode are excluded. Restoring it does not execute the solver binary.

From the repository root:

```sh
python3 research/ic_candidate_tournament_20260915/evidence/restore.py \
  --archive ic-producer-admission-20260925 --out /tmp/ic-producer-replay
python3 research/ic_candidate_tournament_20260915/goal_20260924/producer-admission/audit_transport.py \
  /tmp/ic-producer-replay/goal_20260924/producer-admission/raw/36106499204/ic-producer-both-36106499204-1
```

Repeat the last command with `scaled` and `pairinv` in place of `both`. It requires
Python, not Valgrind or a Rust rebuild. Instruction/native values in the receipts
are fixed-vector diagnostics, not a paired comparison across these CI hosts.

## Remaining campaign gates

Integrate these admission records into the existing development and promotion
drivers, measure instrumentation/admission overhead against original sources,
qualify the strongest compatible IC and matched rho reference, and implement
public-point input plus separately bounded single-target online timing. Native
per-phase wall times and the online interval are still unknown here. Retain the
active goal's complete cold instruction/native acceptance gates; do not substitute
an online-only or amortized result for that objective.

Then seal the five-or-more-cell panel, fresh confirmation targets and familywise
rule before running up to three improvement rounds. Preserve diverse mechanisms
and measured combinations, failures, confirmation and replay. The active goal
remains incomplete until the bounded campaign and its scoreboard/evidence
requirements are fulfilled.
