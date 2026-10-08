# Conditional preregistration: paired n13 O-aware sparse versus dense SAT

Status: **HELD**, frozen source/input gate only; no paired solver outcome.
The independently validated sparse sequential-one-hot/allowed-output CNF
merged in [PR #786](https://github.com/aburan28/crypto/pull/786). This
successor to [PR #787](https://github.com/aburan28/crypto/pull/787) pins its
merged evidence and runner before any timed child. Independent review and
hash-only CI are required before the separately archived measured attempt.
[PR #785](https://github.com/aburan28/crypto/pull/785) is a historical
single-encoding baseline, not the paired reference for this comparison.

## Identical decision problems and model semantics

Use exactly #785's 32 fixed n13-m5 full-point Q+T targets, in `Q0T0...Q7T3`
order, with the same positive final-state literal for each. The dense base is
merged #781's O-aware rational one-hot CNF. The sparse base is merged #786's
Sinz sequential-one-hot plus allowed-output CNF. Before timing, verify both
base/schema SHA-256 hashes, factor primary groups and x values, all 32
literal/point pairs, source hash/merge ancestry, exact clause/variable header,
and #786's independent branch-complete verdict. Reject any source or label
drift. Both files encode *the same* signed-rational-factor point-target
question, including O-prefix branches; no WDSat same-base m! ordering is
permitted for the distinct rotated slots.

Before any child, the runner must prove process-tree enumeration is permitted on
the execution host; an inaccessible monitor fails closed and archives the
refusal. The `preflight` mode records this gate without launching any
solver or producer. The cold exporter and runner use an exact pinned Python 3.12.8 executable and psutil 7.2.2. The three engines are the exact #785 SHA-pinned installed CryptoMiniSat 5.14.7,
Kissat 4.0.4 and CaDiCaL 3.0.1 binaries with identical CLI flags, single
thread and exit/model parser. Run two tiny SAT/UNSAT interface smokes per
engine before this panel; archive failures. For SAT, decode the five chosen
factor x literals from the corresponding schema, independently derive every
complete rational point fibre using #767's bit-serial/Fermat group law, and
try at most 2^5 sign assignments to exhibit an exact Q+T sum. Do not use
#767's saved witness as a solver certificate. For UNSAT, compare the status
with #767's exhaustive point oracle; there is **no separately checked
DRAT/LRAT proof certificate**. Any timeout, missing status, invalid model or
oracle mismatch is censored/failed, never UNSAT.

## Frozen paired order, caps and complete charging

First run four cold producer children in exact ABBA order: dense-1,
sparse-1, sparse-2, dense-2. Each producer emits the four toy panels
because the merged CLI has no n13-only option. Check the fresh n13
base/schema bytes against the merged archives before querying. Charge the
**full four-panel** producer wall to n13 for each representation; this is a
conservative equal **panel scope**, though dense and sparse perform different work. Keep the two opposite-order
pairs separately; do not choose a favorable export measurement post hoc.
Each exporter has a 180 s wall budget and a 512 MiB sampled process-tree
RSS threshold. The 20 ms sampler kills the process group when it observes an
overrun; an unobserved memory spike or up to one polling interval of wall
overrun is possible, so this is not an OS-enforced memory limit. Preserve stdout/stderr, exact command, source/hash, UTC, wall/CPU/RSS and
failed partial outputs. Semantics audits in #781/#786 remain prior certified
source gates; fresh byte identity makes their verdict applicable to the
new base. A separate audit-inclusive cost may add independent verifier
wall, but the operational portfolio charges exact-point model checks.

For each of the 32 targets, construct the dense and sparse derived CNFs
from their own hash-pinned base by replacing the header clause count with
count+1 and appending the **same** positive target unit. Record full-file
hash/bytes and construction wall for each; never commit 192 base copies.
For each target, rotate engine order by target ordinal modulo 3 as #785
did. Within each `(target,engine)` pair, run adjacent **cold processes**
for the two representations. Dense runs first on even target ordinals;
sparse runs first on odd ordinals. Thus each engine has 16 dense-first and
16 sparse-first target pairs, with no adaptive ordering by target class,
prior solve status or expected hardness. Each cold child reads the entire
appropriate derived input. Process wall cap is 15.0 s, sampled child-tree
RSS cap 2,147,483,648 bytes, poll interval 20 ms; kill the process group
on cap. Archive all 192 raw outputs, exact commands, exit/censor reasons,
CPU, RSS, parse/point-check wall, UTC and source/input digests. OS page cache
is not forcibly flushed; balanced adjacent order limits but does not erase
cache and thermal effects. Source/preflight, both base loads, all query-file
constructions, children and verification are measured rather than free.
The nominal sum of child wall budgets is 4×180 + 192×15 = 3,600 s,
excluding polling/teardown, preflight, setup and verification. A separate
3,900 s parent-wall alarm covers preflight through the panel result, kills
an active child process group on interruption, and retains a failed/censored
receipt and partial raw outputs. Archive replay happens afterward and has
its own CI budget.

For each engine and representation report (1) solver-child wall/CPU alone,
(2) query construction plus child plus point verification, and (3) two
fully charged 32-branch portfolio views, one using exporter pair 1 and
one pair 2, each adding source/preflight and one full producer. Also
report maximum single-child wall separately, all four negative Q full-
branch costs, and the fixed-order first independently verified positive
witness costs for Q0–Q3. The separate audit-inclusive view adds #786's archived independent dense/sparse replay wall to each representation and
is labelled as an earlier, different timed block rather than an operational
per-query requirement. The actual 192-child campaign wall is separate
from the per-engine operational views; no metric may substitute a max
query for a full portfolio. Keep per-Q paired differences and descriptive
intervals separate from the source/export cost.

The correctness gate requires all 192 children to finish with the same
oracle-consistent branch statuses, all SAT witnesses independently replayed,
zero malformed/timeout/over-RSS outcomes, and a passing raw archive replay.
For each engine, call sparse a *robust toy wall improvement* only when its
fully charged sparse/dense ratio is **below 0.90 in both opposite exporter-
order pair views**. If an order pair reverses the ranking, call that engine
inconclusive. An overall sparse next-rung preference requires all three
engines to pass the same gate; otherwise report mixed/inconclusive while
retaining any engine-specific result. This is a bounded engineering
observation, not a 95% population speedup claim. n131 exporter width,
relation rank/yield, full ECDLP S, matched rho and a Certicom logarithm
remain unset.

The successor draft PR pins merged #785/#786/#787 commit ancestry,
source/artifact hashes, binaries, targets, exact commands and caps in
`FROZEN.json` and passes hash-only CI before any smoke or producer. It leaves
the measured decision HELD. After independent review, coordinate a quiet
local host window, retain every failed attempt under a unique path, and
commit raw archive, independent replay, analysis, scoreboard and decision
in a linked outcome PR. Exact-head CI/review and guarded merge remain the
final evidence gate.
