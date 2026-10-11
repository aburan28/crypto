---
title: "sect113r1: additive three-node GHS/Hess structural screen"
author: "crypto research harness"
date: "2026-10-07"
lang: en-US
toc: true
toc-depth: 2
geometry: margin=0.72in
fontsize: 10pt
mainfont: "DejaVu Sans"
monofont: "DejaVu Sans Mono"
colorlinks: true
linkcolor: blue
urlcolor: blue
---

# Status and question

**ADDITIVE PRE-ADMISSION STRUCTURAL DIAGNOSTIC.** This addendum asks whether
the direct classical GHS/Hess structural screen identifies a low-genus descent
candidate on `sect113r1`
(`icv1-f2m113-tm122610772499221213-97df4ac6`) or either of the two rational
degree-5 neighbors already certified by the parent diagnostic.

It does not modify or promote the frozen
[`pre-admission-certificate.json`](../../diagnostics/pre-admission-certificate.json).
That certificate remains a same-implementation pre-admission record whose GHS
field says `NOT_CERTIFIED_MAGIC_NUMBER_UNCOMPUTED`. The three new raw outputs
and their custody metadata are preserved in this directory under the separate
diagnostic identity `SECT113R1-GHS-3NODE-20261007-DIAG1`.

The bounded result is:

> The native structural screener returned magic number 113, type I, and genus
> `2^112 - 1` on all three tested nodes. It found zero candidates within the
> frozen genus bound 64. This is a negative structural result at exactly three
> nodes, not an exhaustive isogeny-class result and not a DLP experiment.

![The prior certified degree-5 routes feed three separate structural screens.
Every screen returns the same high-genus tuple; independent verification and
all attack-transfer obligations remain open.](ghs-screen-flow.svg){width=98%}

# Frozen objects and method

All curves use the binary Weierstrass model

```text
y^2 + x*y = x^3 + a*x^2 + b
```

over `GF(2^113)` in the polynomial basis
`z^113 + z^9 + 1`, encoded as
`0x20000000000000000000000000201`. The shared coefficient is
`a=0x3088250ca6e7c7fe649ce85820f7`.

The three frozen node records are:

```text
source      icv1-f2m113-tm122610772499221213-97df4ac6
            EC1N113Csect113r1hf529f17bd191
            b=0xe8bee4d3e2260744188be0e9c723
codomain A  icv1-f2m113-tm122610772499221213-fd54d0eb
            EC1N113Crbh921ab2cd913f
            b=0x109267245489e254e8f14002629a1
codomain B  icv1-f2m113-tm122610772499221213-5de3030f
            EC1N113Crbh2ee3581888f2
            b=0x162b1a595685a1387c82647bf44cf
```

The complete UIDs, exact commands, source hashes, and output hashes are in
[`manifest.json`](manifest.json). The endpoint identities and degree-5 route
evidence originate in the frozen parent
[`isogeny-routes.json`](../../isogeny-routes.json); this addendum does not
re-prove those maps.

The release binary enumerated every nontrivial factorization `113=n*l`. Since
113 is prime, the only row is `(n,l)=(113,1)` and the only field tower is

```text
F_2 <= F_(2^1) <= F_(2^113).
```

For each node the implementation computes the full Hess magic number as an
`F_2`-rank of Frobenius-orbit generators, selects its recorded type branch,
and derives the genus and Artin--Schreier cover degree. The genus bound 64 was
frozen on the command line before reading the output.

# Exact observations

| node | magic `m` | type | genus | within bound 64 |
|---|---:|:---:|---:|:---:|
| source (`…97df4ac6`) | 113 | I | `2^112-1` | no |
| codomain A (`…fd54d0eb`) | 113 | I | `2^112-1` | no |
| codomain B (`…5de3030f`) | 113 | I | `2^112-1` | no |

The genus is exactly `2^112-1`; the cover degree is exactly
`2^113=10384593717069655257060992658440192`. The preserved outputs are
[`source.json`](source.json), [`degree5-a.json`](degree5-a.json), and
[`degree5-b.json`](degree5-b.json). A second execution produced byte-identical
JSON for every node.

![All three tested nodes have magic number 113. The dashed line at `m=7`
marks the largest magic number whose type-I genus is at most the frozen bound
64. The editable data are in `ghs-magic-comparison.csv`.](ghs-magic-comparison.svg){width=96%}

The shell timing observation rounded each execution to `0.001` seconds. It is
custody metadata only: there were no repetitions suitable for performance
statistics, no isolated-run protocol, and no baseline/candidate timing
question. No speed claim is made.

# Transfer and feasibility interpretation

The prime subgroup order remains

```text
n = 5192296858534827689835882578830703.
```

For a genus-`g` curve over `F_2`, the Weil bound gives

```text
#Jac(C)(F_2) <= (1 + sqrt(2))^(2g) = (3 + 2*sqrt(2))^g.
```

An exact integer recurrence verifies
`(3+2*sqrt(2))^44 < n`. Therefore a target of genus at most 44 cannot
contain a nonzero image of this prime-order subgroup: a homomorphism from a
prime-order group is either injective or zero. This is a necessary capacity
gate, not a sufficient attack test. The observed genus `2^112-1` is far above
that lower gate and far above the frozen structural bound 64.

The narrow supported conclusion is
`NO_LOW_GENUS_GHS_CANDIDATE_AT_THE_THREE_TESTED_NODES`. It does not establish
any of the following:

- an independent implementation replay of the Frobenius rank or type rule;
- absence of a special representative elsewhere in the isogeny class;
- an explicit smooth descended curve or Jacobian;
- preservation of the order-`n` subgroup, inverse recovery, or a verified DLP;
- an end-to-end cost below the source curve's generic reference; or
- a prevalence estimate for the isogeny class.

The source curve remains class-weak at its already documented legacy generic
baseline of `2^56.325748...` expected group operations. This screen found no
additional GHS weakness on the two tested valid neighbors; it neither changes
nor strengthens the separate conditional unchecked-input finding involving a
singular, non-isogenous cubic.

# Reproduction and validation

The tracked worktree was clean. The exact producer revision was:

```text
1f0d0bd7863c78ac9d7291b84e78ed937e8c4e2f
```

The release binary SHA-256 was:

```text
498b711179ec8a9dd9e7a12528b2025b1b7f7a89f210c28ccc6f625dc899f9ca
```

For each node, replace `<B>` and `<OUT>` below with the values frozen in the
manifest:

```sh
ghs_screen \
  --degree 113 \
  --modulus 0x20000000000000000000000000201 \
  --a 0x003088250CA6E7C7FE649CE85820F7 \
  --b <B> \
  --genus-bound 64 \
  --output <OUT>
```

Focused release tests at that revision reported:

```text
filter                 passed  failed
library ghs_screen          3       0
binary ghs_screen           2       0
library ec_trapdoor         9       0
library sect113r1          11       0
```

The additive integration verifier checks the three raw hashes, reruns the
structural arithmetic through the native API, verifies the exact result tuple,
and proves the genus-44 capacity inequality using integer arithmetic. It is a
focused same-codebase verifier, not source-independent reproduction. The
release command

```sh
cargo test --offline --locked --release --test sect113r1_ghs_evidence -- --nocapture
```

reported 2 passed, 0 failed, 0 ignored. The manifest pins the verifier source,
test binary, toolchain, environment, and complete command.

# Graph and dashboard scope

The canonical index-calculus scoreboard, progress timeline, leaderboard,
curve registry, cover catalog, and prior study-local visual were checked for
applicability. None changes: this addendum adds no curve identity, no isogeny
edge, no admitted IC measurement, no matched rho comparison, and no
performance result. The two visuals in this directory are the additive record
for the new structural trait. The parent report and its frozen diagram remain
unchanged.

# References

1. P. Gaudry, F. Hess, and N. Smart, “Constructive and destructive facets of
   Weil descent on elliptic curves,” *Journal of Cryptology* 15 (2002).
2. F. Hess, “Generalising the GHS attack on the elliptic curve discrete
   logarithm problem,” *LMS Journal of Computation and Mathematics* 7 (2004).
3. A. Menezes and E. Teske, “Cryptographic implications of Hess' generalized
   GHS attack,” *Applicable Algebra in Engineering, Communication and
   Computing* 16 (2006).

These references provide theorem context. The three output files and
`manifest.json` are the evidence for this bounded diagnostic.
