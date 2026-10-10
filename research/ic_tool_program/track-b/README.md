# Track B's implementations, on record before their measurements

**Track B's steps are declared, and their code is written but not yet
measured.** Each step's protocol and frozen cases are merged. Each
step's code reaches `src/` only through its own results pull request,
after its declared measurements run on the programme's timed queue. That
queue is B0, B1, B3, B2, B2b, B7a and B3b, after R03 and R02b, and then
B4.

Until then the code lives on local branches. This directory keeps it on
record, so that losing a container loses no code. AGENTS.md asks for
this: a local branch must not be the only record of the work.

## What is here

- [`stack-20261001.bundle`](stack-20261001.bundle): a git bundle of every
  local Track B branch, against main at `8b5fdd01`.
- `stack-20261001.bundle.sha256`: its SHA-256.
- [`stack-20261001.txt`](stack-20261001.txt): each branch's tip, and every
  commit in the bundle, oldest first.
- [`stack-20261001-b4.bundle`](stack-20261001-b4.bundle), its
  `.sha256` and [`stack-20261001-b4.txt`](stack-20261001-b4.txt): B4's
  branch, on `b3b-int` at `eaa842aa`, which the first bundle holds.

The branches:

| branch | tip | what it holds |
|:--|:--|:--|
| `b0-port` | `216dbcd8` | B0: the refusals and fixes for the survey's defects, with conformance suite v1, on v0′ (`c1a2e5f8`) |
| `b1-local` | `4deec0de` | B1, on B0: schema v2, the checks for binary instances, the Koblitz pipelines on imported instances |
| `b3-int` | `911b24d7` | B3, on B1: two-word binary fields, `rho-koblitz` on them, and `solve: rho` |
| `b2-int` | `474659af` | B2, on B3: prime and extension fields, the other curve forms and importers, `rho-bignum`, estimates |
| `b2b-int` | `5c93e142` | B2b, on B2: `kic` on curves over a subfield, `kic` alone, order certificates |
| `b7a-int` | `58cada36` | B7a, on B2b: the F1 sampled level for `kic` and every rho |
| `b3b-int` | `eaa842aa` | B3b, on B7a: the index calculus on two-word fields at F0, with its amendments 1–2 |
| `b4-int` | `8a2519e6` | B4, on B3b: binary fields of three to nine words, `rho-koblitz` and `kic` on them, and v1's three rules past one word (in `stack-20261001-b4.bundle`) |
| `b0-local` | `a1187008` | an earlier version of B0, which `b0-port` supersedes |
| `b3-local` | `65465f68` | an earlier version of B3, which `b3-int` supersedes |

The first seven form one stack, each on the one before, and `b4-int`
continues it. The two earlier
versions are outside it, and are kept so that nothing is lost.

## Restoring it

From a clone of this repository, in this directory:

```
sha256sum -c stack-20261001.bundle.sha256
git bundle verify stack-20261001.bundle
git fetch stack-20261001.bundle 'refs/heads/*:refs/heads/track-b/*'
sha256sum -c stack-20261001-b4.bundle.sha256
git bundle verify stack-20261001-b4.bundle
git fetch stack-20261001-b4.bundle 'refs/heads/*:refs/heads/track-b/*'
sha256sum -c stack-20261005-main.bundle.sha256
git bundle verify stack-20261005-main.bundle
git fetch stack-20261005-main.bundle 'refs/heads/*:refs/heads/track-b/*'
```

The second bundle builds on the first's `b3b-int`, so it is fetched
second. The third builds on both and on main's `995ea207`, so it is
fetched last.

The bundle's prerequisites are commits on main, so the clone needs
main's history up to `8b5fdd01`. The seven commits the bundle builds
on, all on main, are listed in `stack-20261001.txt`.

## The arms on main (2026-10-05)

The programme re-based on main's head, `995ea207`, in R07. Track B's
measurements therefore run there, on arms rebuilt from the branches
above. [`stack-20261005-main.bundle`](stack-20261005-main.bundle), with
its `.sha256` and [`stack-20261005-main.txt`](stack-20261005-main.txt),
holds them. Its prerequisites are main's `995ea207` and the eight
step tips, which the two bundles above hold.

**Each arm is built on the one before**, as the chain's amendment
requires. Each is one merge of its step's tip into the arm before it,
from main's head:

| arm | commit | what it merges |
|:--|:--|:--|
| `tbarm-B0` | `2177afa8` | `b0-port` into main's `995ea207` |
| `tbarm-B1` | `ac74627e` | `b1-local` into arm B0 |
| `tbarm-B3` | `ff3b9624` | `b3-int` into arm B1 |
| `tbarm-B2` | `603c6317` | `b2-int` into arm B3, and the port below |
| `tbarm-B2b` | `1aed30b0` | `b2b-int` into arm B2 |
| `tbarm-B7a` | `570d4e9a` | `b7a-int` into arm B2b |
| `tbarm-B3b` | `98320f2f` | `b3b-int` into arm B7a, and the port below |
| `tbarm-B4` | `d5761c13` | `b4-int` into arm B3b, and the port below |

**The merges' conflicts are resolved by rule**, by one script, so that
every arm is resolved the same way:
- **Main's `koblitz_wide`.** Main has its own `koblitz_wide.rs`
  (ecbench's wide Koblitz curve, #1360), so main's module keeps the
  name. The stack's two-word curve and matched rho (B3) becomes
  `koblitz_two_word.rs`, unchanged, and every use of it is renamed.
- **A file only the stack changes.** Where the arm's only edits are
  those renames, the step's own version is taken, renamed.
- **Protocols in conflict.** The union of both sides is taken. In each
  case it is main's text, which already carries the stack's amendments.

**Three port commits** carry what main's changes require. Each is the
same change as on the merged stack (`6a787e07`):
- **At B2:** `ic`'s estimates find each size by its ICV1 slug. Main's
  v1 and v2 baselines carry no `size` label. Rho's figure comes from
  the newest baseline that measured it (v0's).
- **At B3b and B4:** the two-word and W-word pair tables take R05's
  presence filter, as main's one-word table does. Sameness with the
  one-word pipeline then holds on main.

**The check.** Arm B4's tree is byte-identical to the stack merged
with main in one step, `52f87c9b` plus `6a787e07`. On that tree:
- the library tests of every Track B module and `ic`'s own tests pass
  (209 and 45);
- the conformance suite at every step through B4 passes 95 of 95 cases,
  natively (`icprog conformance`).

**Each arm, checked at its own level.** Each arm was built with its own
commit embedded (`IC_BUILD_COMMIT`). Its library tests (every Track B
module's, and `ic_boundary`'s), `ic`'s own tests and the conformance
suite through its step all pass. These are checks of the arms, not the
steps' measurements, which run again under each protocol.

| arm | library tests | `ic`'s tests | conformance |
|:--|--:|--:|--:|
| B0 | 202 | 10 | 8 of 8 |
| B1 | 203 | 17 | 31 of 31 |
| B3 | 208 | 17 | 37 of 37 |
| B2 | 211 | 30 | 57 of 57 |
| B2b | 212 | 33 | 69 of 69 |
| B7a | 213 | 40 | 76 of 76 |
| B3b | 223 | 41 | 83 of 83 |
| B4 | 209\* | 45\* | 95 of 95 |

\* On the stack merged in one step, whose tree is arm B4's; its library
run left out `ic_boundary`'s tests.

The reports, the script that ran them, its log and the binaries' SHA-256
are in [`arms-check-20261005/`](arms-check-20261005/), with their
`SHA256SUMS`.

## B5a's arm on main (2026-10-06)

B5a's re-declaration
([`../rounds/B5a-extension-fields/PROTOCOL.md`](../rounds/B5a-extension-fields/PROTOCOL.md))
puts its arm after B4's, built the same way.
- **The arm:** `tbarm-B5a` (`97ade3d0`), one merge of B5a's branch
  (`b5a-int`, `538afdce`) into arm B4, by the rules above. Two conflicts
  were resolved by hand:
  - `src/cryptanalysis/mod.rs`: the arm's `koblitz_two_word` kept beside
    B5a's modules;
  - `src/cryptanalysis/gaudry_cubic.rs`: main had added
    `Fp3::with_cube_nonresidue` for the field B5a's `Fp3::with_c` builds,
    under the same conditions, so `with_c` returns its result.
- **On record** in [`stack-20261006-b5a.bundle`](stack-20261006-b5a.bundle)
  (its SHA-256 beside it), with B5a's two commits and the merge; arm B4,
  its prerequisite, is in `stack-20261005-main.bundle`. Its listing is
  [`stack-20261006-b5a.txt`](stack-20261006-b5a.txt). Restore it as
  above, after the earlier bundles.
- **Checked at its level,** as the other arms were, with its commit
  embedded: 259 library tests (`fpk_curve`, `curve_id`, `gaudry_cubic`,
  `rho_bignum` and the Koblitz pipelines'), 50 of `ic`'s and 8 of
  `tests/curve_id.rs` pass, and the conformance suite through B5a passes
  110 of 110 cases (C050 superseded by C103). The report, script, log
  and binaries' SHA-256 are in
  [`arms-check-20261006-b5a/`](arms-check-20261006-b5a/).

## What this is not

- **Not accepted code.** No step here has been measured. Each is judged
  by its own protocol, and an amendment or a failed measurement can
  still change it.
- **Not a baseline,** and not a claim about speed or correctness.
- **Not a change to `src/`.** Main's tool is unchanged by this
  directory.

When a step's results pull request merges its code into `src/`, the
bundle stays as the record of what was measured. A later bundle goes
in a new file; none is overwritten.
