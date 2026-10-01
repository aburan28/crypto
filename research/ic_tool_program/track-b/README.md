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
```

The second bundle builds on the first's `b3b-int`, so it is fetched
second.

The bundle's prerequisites are commits on main, so the clone needs
main's history up to `8b5fdd01`. The seven commits the bundle builds
on, all on main, are listed in `stack-20261001.txt`.

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
