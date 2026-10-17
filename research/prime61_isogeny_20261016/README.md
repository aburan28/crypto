# Prime61 horizontal-isogeny factor-base panel

2026-10-16. **Partial requested experiment:** completed every outgoing
kernel-distinct rational isogeny of prime degree at most 61 from the source;
the entire horizontal class has not been exhausted. The earlier
[source-only baseline](../prime61_factor_base_20261016/README.md) is separate.
This is a null relation-stage result, not a measured improvement.

| Requested requirement | Status and evidence |
| --- | --- |
| Use the repository code at a requested large field | Complete at the prime field `p=2^61-1`; unchanged native `prime_fast` and `isogeny_walk` backends |
| All horizontal isogenies | Partial: 15 outgoing prime-degree kernel edges through 61, one hop only; full-class coverage remains [blocked](../../problems/prime61-horizontal-class-memory.md) |
| Several vertical isogenies, if they exist | Resolved for separable isogenies defined over this base field: none exist in this class, since `f_pi=f_E=1` |
| Test mapped instances | Complete for the finite panel: 16 workloads, 12 paired streams each; every stream reached its cap without a relation or solution |
| Reproduce the screenshot's exact 1,000/2,000-candidate ranking | Still open: [original-driver problem](../../problems/screenshot-factor-base-driver-unlocated.md); this native panel uses three inverse pairs and does not recreate that search |

The phrase “all horizontal” could mean all class vertices or outgoing edges
under a degree bound. No clarification was received during execution. The
finite pilot does not silently replace either an unbounded map enumeration
or an exhaustive class experiment.

## Frozen inputs and identity

The source library is `de0c9391480155ebdd5a8495d2c79e4508d4795d`.
The protocol was committed as `596cae7d`, with class-count audit
`14298786`, degree-61 extension `da85e554`, and exact class-order certificate
`95ec4bec`, all before their respective runs. The initial transport driver
was committed as `764d53c0`; `95ec4bec` adds the actual measured EC1 records.
It does not change collection counts. See [PROTOCOL.json](PROTOCOL.json)
and the three named addenda in this directory.

```
p = 2305843009213693951
y^2 = x^3 + 2305843009213693948*x + 111
#E(F_p) = r = 2305843010818082053 (prime), h = 1
G = (1, 1796364306280342617)
Q = (992726311114892094, 558260110083040833)
t = -1604388101
```

This is a public synthetic known-answer target: `Q=[539365399]G`.
The collector receives only `G,Q`, not this scalar. The scalar is used
solely for input and transport replay. Neither this curve nor its target
is a transported version of the screenshot's small-field instance.

Exact measured model, base, transported `G,Q`, EC1 record and full UID for
each workload are the `curve` rows in
[factor-base-l61-identified.jsonl](results/factor-base-l61-identified.jsonl).
The source UID is
`urn:ec-record:1:sha256:21ad10f63efdfb9cece7b48145b24c48a2ee16c9cac75581333f3938db42e6c9`.
Walker destination records use their own generators; their EC1 identities
must not substitute for these measured records. Registration of the new
models in ICV1 remains pending. No unregistered slug is promoted into a
leaderboard, accepted candidate or cross-repository curve link.

## Structure and finite coverage

Exact arithmetic gives

```
Delta_pi = t^2 - 4p = -6649310858224389603
abs(Delta_pi) = 3 * 97 * 113 * 35831 * 5643470711
```

The driver reconstructs the product and deterministically verifies every
factor prime with the existing u64 arithmetic. Each exponent is one and
`Delta_pi=1 mod 4`, so it is a fundamental discriminant. In an ordinary
finite-field isogeny class, `Z[pi]` is contained in every endomorphism order
and `f_E` divides `f_pi`. Consequently `f_pi=f_E=1` throughout this class.
Vertical separable base-field edges would change this conductor, so none
exist here. This says nothing about vertical maps over extension fields or
other source classes. Frobenius/inseparable maps are outside the pilot.

The usual split/ramified/inert classification of primes in the maximal
order gives the complete outgoing prime-degree roster below. All retained
kernels are distinct; their 15 destination j-invariants are also distinct.
These are quotient maps up to codomain field isomorphism, not every
postcomposition with an automorphism or every composite-degree map.

| Prime degree | Horizontal kernel edges | Destination streams |
| --- | ---: | ---: |
| 3 (ramified) | 1 | 12 |
| 7 | 2 | 24 |
| 11 | 2 | 24 |
| 13 | 2 | 24 |
| 19 | 2 | 24 |
| 23 | 2 | 24 |
| 53 | 2 | 24 |
| 61 | 2 | 24 |
| 2 | 0 | 0 |
| 5, 17, 29, 31, 37, 41, 43, 47, 59 (inert) | 0 | 0 |
| Source identity workload | — | 12 |

The degree-2 native rational-kernel solver returned zero. Initial
[class-initial.json](results/class-initial.json) contains the odd-prime
splitting audit; its probable-prime annotations are superseded for the
discriminant factors by the deterministic manifest checks. The degree-13
pilot actually tried the inert degree 5 as a zero-neighbour control.
Other inert odd degrees were excluded from construction by the exact
splitting test, rather than represented as failed kernel searches.

[one-hop-l61/isogeny_routes.json](one-hop-l61/isogeny_routes.json) contains
the exact source and destination models, edge IDs and kernel polynomials.
The separate native verifier passed 16 curves and 15 edges with two audit
points per edge: [receipt](results/walk-l61-verify.json), route SHA-256
`02ac92324b7161ee9516504b2c08db78c961c8dc6573b5e5bc4376aca127ecd3`.
There were no retained edge failures, repeated destinations or construction
cap hits. The earlier degree-13 pilot has its own raw routes and receipt.

## Transport and measured result

[prime61_isogeny_experiment.rs](../../examples/prime61_isogeny_experiment.rs)
rechecks each kernel using the native division-polynomial, squarefree and
cyclic-kernel tests. It evaluates the standard Vélu point formula using
traces in `F_p[T]/h(T)` and applies the recorded codomain isomorphism. This
supports Frobenius-stable kernels without assuming rational generators.
The underlying formula is already implemented in the repository's Vélu
code; this experiment adds quotient-algebra point evaluation, not a new
general identity.

Both `G` and `Q` pass through the same map. Each destination is checked on
curve, in the order-r subgroup, and against the frozen target relation.
There are 64 scalar transport checks per map (960 total), and all 36 signed
pair sums of the transported three-pair source base are compared as a
multiset on each map (540 total). Those transported-base checks are a
control. The measured collector instead builds a fresh native low-x
three-pair base on each codomain, matching the prior source-only protocol.

All 16 workloads use the same 12 seeds, `2026101600..2026101611`, with
20 required relations and 100,000 target trials per stream. From the
native [summary](results/summary.json):

| Recorded quantity | Result |
| --- | ---: |
| Streams / capped streams | 192 / 192 |
| Target trials | 19,200,000 |
| Collected relations | 0 |
| Verified solves | 0 |
| Full operation cost, S, matched-rho cost, speedup | Unknown |

The driver independently re-adds any returned relations and checks any
recovered scalar; there were no large-field relations to replay. A second
identified execution is byte-identical to the first. This is same-host
reproducibility, not independent-host validation. The original degree-61
run's collection rows also match byte for byte after adding measured EC1
metadata. The 96 degree-13 pilot streams and repeated degree-61 streams
are historical controls and are not pooled into the 192-stream headline.
The earlier 20-bit positive control solved and replayed all its relations;
see the separate baseline report, not a new control claimed here.

These counters describe a capped relation stage. No full pipeline solved,
and construction, transport, setup, complete primitive accounting and
strong-rho references have not been priced. Raw walker files include
unisolated diagnostic timings; they are not admitted performance figures.
Canonical performance dashboards and leaderboard graphs stay unchanged
because no complete-cost result or accepted candidate is added.

## Full-class lower bound and blocker

PARI/GP 2.15.2, a compiled native arithmetic backend, independently returned
`545296208` from `qfbclassno`, `quadclassunit`, and the BNF class-group
calculation, with conditional group structure `[68162026,2,2,2]`.
**The full class number is not certified:** `bnfcertify` hit the frozen
300-second process cap (exit 124). The Shanks routine alone is also not
reliable as a proof at this discriminant size. Preserve this number only
as agreeing, GRH-conditional class-group output.

A separate exact certificate avoids that assumption for a lower bound.
Let `A=Qfb(7,5,237475387793728201)` and `n=68162026=2*34081013`.
Native exact reduced-form composition verifies `A^n` is principal, while
`A^(n/2)` and `A^(n/34081013)` are not. Both prime factors are verified.
Thus A has exactly order n: any proper divisor of n omits at least one
prime factor to its full multiplicity, contradicting the corresponding
nonprincipal power. Hence the class group has at least 68,162,026 elements.
The ordinary maximal-order class-group action gives at least this many
base-field isomorphism classes in the fixed trace class. This is a
standard CM correspondence, not an enumeration or a certified full count.

The existing walker retains `Vec<Node>` and `Vec<Edge>`. Native `sizeof`
is 240 bytes per node and 152 per edge. Even a minimum spanning record
of n vertices and n−1 edges requires **26,719,514,040 bytes in headers
alone**, exceeding this container's **17,179,869,184-byte (16 GiB) limit**.
This excludes every owned allocation, kernel coefficient and hash table.
The first preregistered node-only estimate was 16,358,886,240 bytes and
did not itself exceed the limit; the subsequent edge-header refinement is
recorded separately rather than rewriting that evidence. It establishes a
resource blocker for this in-memory backend, not mathematical impossibility
or a bound on an external-memory implementation.

See [class_order_certificate.gp](class_order_certificate.gp), its
[successful exact transcript](results/class-order-certificate.txt), and
[native memory record](results/min-class-node-edge-memory.json).
The first interactive GP attempt used reserved variable `I`, failed and
nevertheless printed symbolic output and a misleading PASS. Its stdout
and error transcript are retained as `class-order-initial-invalid.*` and
are **invalid evidence**. The replacement uses local variables, type and
value checks, and `iferr(...quit(1))`; its stderr is empty and exit is zero.

The full class remains open in the [named problem](../../problems/prime61-horizontal-class-memory.md)
and [blocked backlog card](../../backlog/XP-20261016-prime61-horizontal-class/card.json).
No full traversal was launched, no remote fleet was scheduled, and no
full-class factor-base result is claimed.

## Verification, proof status and reproduction

- New example unit test: quotient-trace degree-3 evaluation matches the
  existing independently implemented small-degree map, including infinity
  and the kernel pole. It passes. Native walk verification and all recorded
  transport checks pass. The native summarizer checks the edge histogram,
  stream roster, seeds, caps and counters before emitting the summary.
- The unchanged source library's earlier release suite passed 4,154 tests
  with 117 ignored; that existing log is in the baseline study. This is not
  reported as a fresh whole-library run after the panel.
- Strict Clippy first failed on a pre-existing mismatched-width expression
  in `src/symmetric/modes/pmac.rs`. The scoped rerun permits that one lint;
  it passed, as recorded in the retained log. No unrelated library
  source is changed to hide it. Failed example compilation attempts and
  the initial missing-Clippy-component error are retained too.
- Lean is not installed in this environment. No checked Lean theorem is
  claimed. The native squarefree/primality and reduced-form certificates
  are arithmetic evidence, not a Lean formalization of CM classification,
  the class-group action, or Vélu's theorem. Those structural bridges
  remain unformalized here. This report invokes those established theorems
  rather than proposing a new polynomial identity; the polynomial-only
  `identity.certificate/v1` schema does not certify a quadratic-form class
  order or the CM correspondence.
- `conductor check` was attempted before edits and again for the report
  scope; the command is unavailable. The workflow file is not tracked.

[reproduce.sh](reproduce.sh) is thin shell orchestration of native Rust
builds, construction, route verification, capped collection, replay and
summary validation into a new directory. GP scripts reproduce the separate
class audit and exact order certificate; use an outer 300-second timeout.
The native `summary.rs` and `resources.rs` sources are included. Source,
binary, route and evidence SHA-256 lists are in `results/`; the older
`source-and-input.sha256` intentionally describes the initial driver, not
the later EC1 metadata addition. Package hashes pin the native PARI build.

The complete reproduction script exited zero; its fresh records and native
summary are retained under `results/script-reproduction/`, with the build
and test log beside them. The class-order transcript also replayed exactly.
The final resource reader queries `/sys/fs/cgroup/memory.max`; its output
matches the earlier header-bound record. All three PDF pages were rendered
and visually inspected after correcting a Typst path/comment interaction.

Host: Linux x86_64, AMD EPYC 9V74, five visible CPUs, 16 GiB container
limit. The current host capsule is [host.json](results/host.json).
No isolated wall measurement or performance ratio is inferred.
Editable visual and PDF: [coverage.svg](reports/coverage.svg),
[report.typ](reports/report.typ), [report.pdf](reports/report.pdf).
The graph is generated by `visual.rs` from the validated summary; compile
it with the same native `serde_json` dependency linkage as `summary.rs` in
`reproduce.sh`, then redirect its output to `reports/coverage.svg`.
Typst compiles `reports/report.typ` to the adjacent PDF. Python was used
only as the document-renderer binding, never for research computation.

References: [Sutherland, Isogeny volcanoes](https://arxiv.org/abs/1208.5370)
for ordinary CM orders, rational volcano levels and class actions;
[PARI arithmetic manual](https://pari.math.u-bordeaux.fr/dochtml/html-stable/Arithmetic_functions.html)
for `qfbclassno`, `quadclassunit`, and exact binary quadratic forms;
[PARI number-field manual](https://pari.math.u-bordeaux.fr/dochtml/html-stable/General_number_fields.html)
for conditional BNF output and `bnfcertify`.
