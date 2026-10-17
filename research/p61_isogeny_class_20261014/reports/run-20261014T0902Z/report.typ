#let navy = rgb("#17324d")
#let blue = rgb("#2367a8")
#let green = rgb("#15803d")
#let amber = rgb("#b45309")
#let slate = rgb("#475569")
#let pale = rgb("#f3f6fa")
#let rulecolor = rgb("#cbd5e1")

#set document(title: "p61 selected-class isogeny census", author: "Cryptanalysis research pipeline")
#set page(
  paper: "us-letter",
  margin: (x: 0.68in, top: 0.66in, bottom: 0.64in),
  header: context {
    if counter(page).get().first() > 1 {
      set text(size: 7.5pt, fill: slate)
      grid(columns: (1fr, auto), [p61 selected-class census], [Evidence report - 2026-10-14])
      line(length: 100%, stroke: 0.5pt + rulecolor)
    }
  },
  footer: context {
    set text(size: 7.5pt, fill: slate)
    line(length: 100%, stroke: 0.5pt + rulecolor)
    align(right)[Page #counter(page).display()]
  },
)
#set text(font: "Avenir Next", size: 9.2pt, fill: rgb("#172033"))
#set par(leading: 0.63em, justify: true)
#set heading(numbering: "1.")
#show heading.where(level: 1): it => block(above: 13pt, below: 7pt)[#set text(size: 15pt, weight: "bold"); #it]
#show heading.where(level: 2): it => block(above: 10pt, below: 5pt)[#set text(size: 11.5pt, weight: "bold"); #it]
#show link: set text(fill: blue)

#let tag(label, color) = box(
  inset: (x: 7pt, y: 3pt),
  radius: 3pt,
  fill: color.lighten(88%),
  stroke: 0.7pt + color,
  text(size: 7.6pt, weight: "bold", fill: color)[#label],
)

#let callout(title, body, color: blue) = block(
  width: 100%,
  inset: 10pt,
  radius: 5pt,
  fill: color.lighten(93%),
  stroke: (left: 3pt + color),
)[
  #text(weight: "bold", fill: color)[#title]
  #v(3pt)
  #body
]

#align(center)[
  #v(0.25in)
  #text(size: 25pt, weight: "bold", fill: navy)[p61 selected-class isogeny census]
  #v(7pt)
  #text(size: 14pt, fill: blue)[Certified launch evidence and completion boundary]
  #v(19pt)
  #rect(width: 84%, height: 1.4pt, fill: blue)
  #v(17pt)
  #grid(
    columns: (auto, 1fr),
    column-gutter: 14pt,
    row-gutter: 7pt,
    align: (right, left),
    [#text(weight: "bold", fill: slate)[Field]], [`p = 2^61 - 1 = 2305843009213693951`],
    [#text(weight: "bold", fill: slate)[Selected root]], [`p61-i0025`, trace `-704859453`],
    [#text(weight: "bold", fill: slate)[Class order]], [`358094304`, factors `[179047152, 2]`],
    [#text(weight: "bold", fill: slate)[Run]], [`p61-class-pilot-v2-20261014T0902Z`],
    [#text(weight: "bold", fill: slate)[Evidence date]], [2026-10-14],
  )
  #v(21pt)
  #tag("CLASS ORDER CERTIFIED", green)
  #h(8pt)
  #tag("PREFIX AUDIT PASSED", amber)
  #h(8pt)
  #tag("FULL CENSUS INCOMPLETE", slate)
]

#v(24pt)
#callout(
  "Scope",
  [One selected ordinary isogeny class over the field. This report does not claim enumeration of every curve or every isogeny class over the field.],
  color: navy,
)

#v(12pt)
#align(center)[#text(size: 8pt, fill: slate)[Source commits: producer `c9bee1b2`; auditor `1e9f40817`; parallel replay candidate `d06472954`; certification runner `d2f5f49a`.]]

#pagebreak()

= Decision

PARI `bnfcertify(bnf)` with flag 0 returned `1` after 10,365.337 seconds. The
dedicated receipt unconditionally certifies the selected fundamental
discriminant `-8726545188371316595` with class number `358094304` and cyclic
factors `[179047152, 2]`.

The complete bounded pipeline has been exercised on an immutable 300,672-record
prefix: generation, independent mathematical replay, checkpoint and hash-chain
validation, dual-evidence ledger sealing, external sort, and uniqueness audit.
The certified audit found 300,672 unique j-invariants and zero duplicate records.

#callout(
  "Current boundary",
  [Coverage is 300,672 of 358,094,304 vertices, approximately 0.0839645%. The receipt status is `partial_census_pass`; `complete_class_coverage` is false. The degree 5/13/37/59 adjacency pass is separate and pending.],
  color: amber,
)

#v(7pt)
#figure(
  image("evidence-pipeline.svg", width: 100%),
  caption: [Certified inputs, completed prefix gates, and the remaining full-census and adjacency gates.],
)

#pagebreak()

= Requested outcome and evidence state

#set table(inset: 5pt, stroke: 0.45pt + rulecolor)
#table(
  columns: (1.65fr, 1.05fr, 2.6fr),
  table.header(
    table.cell(fill: navy)[#text(fill: white, weight: "bold")[Component]],
    table.cell(fill: navy)[#text(fill: white, weight: "bold")[State]],
    table.cell(fill: navy)[#text(fill: white, weight: "bold")[Evidence boundary]],
  ),
  [Exact selected-class order], [#tag("CERTIFIED", green)], [Flag-0 certificate returned `1`; receipt SHA-256 `260382fb...f9643c9`.],
  [Bounded-memory producer], [#tag("IMPLEMENTED", green)], [Committed Rust producer with checkpoints, chained records, compressed segments, and failure retention.],
  [Independent replay], [#tag("PASS", green)], [All four segments mathematically replayed with distinct order-audit settings.],
  [External uniqueness], [#tag("PARTIAL PASS", amber)], [300,672 records sorted; 300,672 unique invariants; zero duplicates; incomplete class range.],
  [Full vertex census], [#tag("INCOMPLETE", slate)], [300,672 of 358,094,304 vertices.],
  [Selected-degree adjacency], [#tag("PENDING", slate)], [2,506,660,128 expected directed degree 5/13/37/59 records after the full join.],
  [All field-wide curves/classes], [#tag("NOT ENUMERATED", slate)], [A different and much larger objective than completing this selected ordinary class.],
  [M7 large-field use], [#tag("NOT ADMISSIBLE YET", amber)], [Requires complete graph evidence, frozen labels and masks, and leakage-safe successor controls.],
)

== Mathematical identity

#table(
  columns: (1.45fr, 2.55fr),
  table.header(table.cell(fill: pale)[#text(weight: "bold")[Quantity]], table.cell(fill: pale)[#text(weight: "bold")[Value]]),
  [Prime], [`2305843009213693951`],
  [Trace], [`-704859453`],
  [Fundamental discriminant], [`-8726545188371316595`],
  [Conductor], [`f_pi = 1`],
  [Class number], [`358094304`],
  [Cyclic factors], [`[179047152, 2]`],
  [Certification], [`bnfinit(polynomial, 1)` then `bnfcertify(bnf)`],
  [PARI flag/result], [`0 / 1`],
  [PARI version], [`2.17.2`],
  [Wall time], [`10365.721123667 s`],
)

#pagebreak()

= Certified prefix generation and replay

The traversal has 18 lanes over 19,894,128 columns. The current prefix contains
two 32-column pilot segments, one 256-column expansion, and one 16,384-column
certified expansion.

#table(
  columns: (1.2fr, 0.72fr, 0.9fr, 0.95fr, 0.95fr),
  align: (left, right, right, right, right),
  table.header(
    table.cell(fill: navy)[#text(fill: white, weight: "bold")[Segment]],
    table.cell(fill: navy)[#text(fill: white, weight: "bold")[Vertices]],
    table.cell(fill: navy)[#text(fill: white, weight: "bold")[Generate]],
    table.cell(fill: navy)[#text(fill: white, weight: "bold")[c9 replay]],
    table.cell(fill: navy)[#text(fill: white, weight: "bold")[d064 replay]],
  ),
  [`x = 0..32`], [576], [13.415 s], [14.778 s], [23.204 s],
  [`x = 32..64`], [576], [13.421 s], [14.672 s], [22.693 s],
  [`x = 64..320`], [4,608], [21.266 s], [46.939 s], [33.460 s],
  [`x = 320..16,704`], [294,912], [858.472 s], [not run], [735.772 s],
  [#text(weight: "bold")[Total]], [#text(weight: "bold")[300,672]], [#text(weight: "bold")[906.574 s]], [#text(weight: "bold")[76.389 s†]], [#text(weight: "bold")[815.129 s]],
)

#v(6pt)
The c9 total excludes the largest segment. The committed `d06472954` parallel
replay candidate passed all four immutable segments and reproduced their
record-chain, checkpoint, and segment digests.
It recorded 23.204 s and 22.693 s on the short segments: 1.570x and 1.547x the
earlier replay times. On the 4,608-record expansion, the first comparison was
33.460 s for d064 versus 46.939 s for c9, a 1.403x c9/d064 ratio. A second c9
control took 215.940 s during several release compiles. Every identity agrees,
but this variance precludes a general speedup claim; controlled isolated ABBA
runs remain required.

On the 294,912-record expansion, generation took 858.472 s and d064 replay took
735.772 s (343.531 and 400.820 records/s). No matched c9 replay was run there,
so this establishes sustained operation and correctness, not a speedup.

#v(8pt)
#figure(
  image("pilot-timings.svg", width: 72%),
  caption: [Startup-inclusive prefix measurements. One observation per bar on a shared, unisolated host.],
)

Extrapolation from the 16,384-column expansion gives about 12.1 days for
generation and 10.3 days for d064 replay, or 22.4 days for those phases alone.
These are diagnostic projections, not schedule commitments; audit,
interruptions, and adjacency are additional.

#pagebreak()

= Certified ledger and durability

The dedicated certificate binds the exact provisional preflight receipt by
SHA-256. Existing segment receipts retain their original preflight binding;
the v2 ledger adds the certificate as a second immutable dependency instead of
rewriting historical evidence.

#table(
  columns: (1.7fr, 1.25fr),
  align: (left, right),
  table.header(table.cell(fill: pale)[#text(weight: "bold")[Audit quantity]], table.cell(fill: pale)[#text(weight: "bold")[Result]]),
  [Ledger range], [`x = 0..16,704`],
  [Records], [300,672],
  [Unique j-invariants], [300,672],
  [Duplicate records/groups], [`0 / 0`],
  [Complete class coverage], [`false`],
  [Status], [`partial_census_pass`],
  [Sorted-index bytes], [4,810,768],
  [Sorted-index SHA-256], [#text(font: "Menlo", size: 6.4pt)[055fb9e007b0b7a261f4850b78450d1d8f306ec56e13c316cebf20514d1793ad]],
)

== Durable storage semantics

The active UGREEN SMB export is mounted at `/Volumes/data`; the requested alias
`/Volumes/UgreenNAS` is not mounted. Artifacts live under
`/Volumes/data/ec-hardness-atlas-20261012/`.

The SMB mount reports `sync_all_unsupported_by_filesystem` for derived ledger,
chunk, index, and receipt files. That limitation is recorded rather than
relabeled. A host `sync` barrier completed, then all terminal digests were
replayed. The certification runner's receipt reports `fsync_completed`.

== Preserved failures and controls

- The first class-number attempt failed during `bnfinit` with PARI's
  low-precision `small_norm` error. Its event log and failure receipt remain.
- The successful runner retried at 256 digits and retains the initial error.
- The flag-1 quotient-only control was interrupted after it remained
  non-decisive; its control receipt remains content-addressed.
- Interrupted segment and SMB per-file-sync attempts remain under
  `p61-class-streaming-v1/incomplete/` and are not admitted.
- The first `d06472954` release build failed during archive creation because
  the system data volume had only 187 MiB free. It remains classified as an
  infrastructure failure in `build-attempts.json`; the corrected build uses a
  scoped target on SSD990-2. After that build passed, only the generated
  `/private/tmp` target was deleted, restoring local headroom to 5.5 GiB.
- The auditor rejects an unconditional label without the dedicated flag-0
  certificate. The false declaration exited nonzero and produced no accepted
  ledger.
- Two default-parallel whole-library release-test attempts remained non-green in the
  pre-existing `jv_cover` module (2 and 3 failures respectively); every exposed
  failure passes in exact isolation and the p61 auditor gate passes 6/6. A full
  low-concurrency run with `--test-threads=1` passed 4,079/4,079, with 114
  ignored, in 2,765.68 seconds. The varying default-parallel failures and
  explicit internal timeout remain in `test-results.json`; those attempts are
  not relabeled as passes. The library is green under the recorded serial mode,
  while default-parallel execution remains non-green on this shared host.

#pagebreak()

= Storage and scale boundary

At the post-expansion observation, the UGREEN export reported 7.3 TiB total, 198 GiB used,
and 7.1 TiB available.

#table(
  columns: (2.1fr, 1fr, 2.1fr),
  table.header(
    table.cell(fill: navy)[#text(fill: white, weight: "bold")[Quantity]],
    table.cell(fill: navy)[#text(fill: white, weight: "bold")[Estimate]],
    table.cell(fill: navy)[#text(fill: white, weight: "bold")[Interpretation]],
  ),
  [All selected-degree directed records], [2,506,660,128], [Seven directed degree 5/13/37/59 neighbors per class vertex.],
  [Raw JSON projection], [about 5.14 TB], [Dominated by the later adjacency representation.],
  [Projection plus 25% headroom], [about 6.43 TB], [Fits only narrowly on the current active volume.],
  [UGREEN available], [7.1 TiB], [Insufficient for redundancy, multiple large fields, and broad intermediate products.],
  [Requested petabyte tier], [not provisioned], [Requires separate object or cold storage plus lifecycle and replication policy.],
)

#callout(
  "Capacity decision",
  [The UGREEN volume can serve as the active staging tier for one selected-class campaign. It is not a petabyte store. A larger program needs packed, indexed graph bundles and separately provisioned durable storage before data generation expands to many large fields.],
  color: amber,
)

Even a provisioned petabyte would not support a literal field-wide census.
Storing only one packed 64-bit j-invariant for each of the `2^61 - 1` possible
j-values already requires `18,446,744,073,709,551,608` bytes (about 18.45 EB,
or 16 EiB), before twists, models, edges, labels, indexes, replicas, or work
space. A petabyte program can cover a very broad selected-class panel, not
every curve over this field.

The current compressed vertex-prefix sample is much smaller than the raw
adjacency projection, but its byte rate must not be transferred to the later
edge corpus because the schemas and record multiplicities differ.

#pagebreak()

= Relationship to M7

The companion M7 implementation uses class-held-out and whole-field-held-out
splits and explicit masks for exhausted, skipped, capped, unverified, and
incomplete-walk observations. The frozen external-field run improves the
ceiling/size baseline, but its preregistered specificity gates fail.

#table(
  columns: (2.3fr, 1fr, 1.7fr),
  table.header(table.cell(fill: pale)[#text(weight: "bold")[Probe]], table.cell(fill: pale)[#text(weight: "bold")[Result]], table.cell(fill: pale)[#text(weight: "bold")[Meaning]]),
  [Primary geometry-residual MLP], [12.8088%], [Baseline test SSE removed; interval 11.0459% to 13.9280%.],
  [Matched geometry permutation], [19.8997%], [A destructive correspondence control improves more than the primary model.],
  [Typed-edge rewiring], [12.5639%], [Topology destruction retains nearly the primary effect.],
  [Node-feature shuffle], [6.1334%], [Node correspondence destruction retains prediction.],
  [Excluded `log2_q` ridge reconstruction], [`R2 = 0.860425`], [Numeric size leakage remains recoverable.],
)

The defensible result is stable predictive correlation for the measured labels,
not isolated geometry-specific incremental information. The p61 graph is an
infrastructure and scale-transfer input; it does not make the M7 specificity
and leakage gates pass.

#pagebreak()

= Completion contract

Full selected-class completion requires all seven stages:

1. Generate all 19,894,128 traversal columns across 18 lanes.
2. Independently replay every admitted segment with distinct audit settings.
3. Seal a ledger covering exactly 358,094,304 records.
4. Externally sort every invariant and prove global uniqueness and exact count.
5. Build and verify the complete degree 5/13/37/59 adjacency join.
6. Attach frozen hardness contexts, censor masks, and provenance.
7. Train and evaluate a leakage-safe M7 successor on held-out classes and
   fields with frozen destructive controls.

#callout(
  "Non-substitution rule",
  [A certified class order does not certify graph coverage. Prefix uniqueness does not certify global uniqueness. A full vertex census does not include adjacency until the join is run. A large graph does not by itself establish a geometry-specific hardness effect.],
  color: navy,
)

Enumerating every curve or every isogeny class over the field is a different,
far larger objective and is not implied by completing this class.

The fail-closed continuation runner has a tested plan beginning at
`x = 16,704`; future chunks are not represented as generated, replayed, or
admitted.

= Reproduction map

- Producer protocol: `docs/p61-isogeny-class-streaming.md`
- Auditor protocol: `docs/p61-isogeny-audit.md`
- Producer CLI: `src/bin/p61_isogeny_class.rs`
- Traversal implementation: `src/cryptanalysis/isogeny_walk/p61_class.rs`
- Auditor CLI: `src/bin/p61_isogeny_audit.rs`
- Fail-closed continuation runner: `research/p61_isogeny_class_20261014/run_remaining_census.py`
- Raw timing record: `timings.json`
- Timing renderer: `plot_timings.py`
- Evidence diagram renderer: `plot_evidence_pipeline.py`
- Release build attempts: `build-attempts.json`
- Release test results: `test-results.json`
- Parallel-verifier release manifest: `parallel-verifier-release-manifest.json`
- Terminal certificate manifest: `certification-terminal-manifest.json`
- Expanded prefix manifest: `certified-prefix-x000000000000-x000000016704-manifest.json`

#pagebreak()

= Digest appendix

#set text(size: 7.4pt)
#set table(inset: 4pt)
#table(
  columns: (1.35fr, 2.65fr),
  table.header(table.cell(fill: navy)[#text(fill: white, weight: "bold")[Artifact]], table.cell(fill: navy)[#text(fill: white, weight: "bold")[SHA-256]]),
  [Provisional preflight receipt], [#text(font: "Menlo", size: 6.2pt)[44e202b8eb9cdad3ef09fb50e416cc8bdcf0942a6f140e24db587a3c4d44493f]],
  [Certification receipt], [#text(font: "Menlo", size: 6.2pt)[260382fbce4805ce84ce00fc77baf27c95243220ae88ca7725e31e4f4f9643c9]],
  [Terminal manifest], [#text(font: "Menlo", size: 6.2pt)[372a304bd30be4aa638264edaf84711d374118637343abd758f0cf00b951b0b1]],
  [Generation producer binary], [#text(font: "Menlo", size: 6.2pt)[47fa827fa9877a91f8b6d496f688f1df6fb18fd25a44b303392ba555be79043f]],
  [Pilot segment `0..32`], [#text(font: "Menlo", size: 6.2pt)[2a29068d66cee5914a5960bff4a13d55dd97ec80d8b9604e4f82c1a6c9c7ddf7]],
  [Pilot segment `32..64`], [#text(font: "Menlo", size: 6.2pt)[70ef588e21200e3c059695ffeae98b29faf0be702b8a01b88b0e12f888bb9fc2]],
  [Pilot completion manifest], [#text(font: "Menlo", size: 6.2pt)[84e597f9f07e674c732e39cabd88d47bf8e902cbeb083492a1e7757bcbf06b8c]],
  [Auditor binary], [#text(font: "Menlo", size: 6.2pt)[bd99626fb402f0e71587b350912f57df021ee5449791ea91bee0a8d2d3917c05]],
  [Parallel replay binary], [#text(font: "Menlo", size: 6.2pt)[81c1a5b8897bcc39fd20c6427e17eb73d9160d838bb001f82810de4272579e82]],
  [Parallel release manifest], [#text(font: "Menlo", size: 6.2pt)[a875ea679bcf04af623b28790baf8f61f0df1f64da7035506bd51da750562268]],
  [Parallel receipt `0..32`], [#text(font: "Menlo", size: 6.2pt)[b24697c9e8fbb53e594864406c3f4c98a6aeb5cd2d349b9c3b3951e0e06f79a3]],
  [Parallel receipt `32..64`], [#text(font: "Menlo", size: 6.2pt)[696d05717de484db6fe2dd92abbf6cf5fd2563e565555d7f560351cb7f5a0a8d]],
  [Certified segment `64..320`], [#text(font: "Menlo", size: 6.2pt)[b0e0180043a3b5aff0a22fa9c0a6151cfbaa0a331b62253f9aa49a97161d8852]],
  [Generation receipt `64..320`], [#text(font: "Menlo", size: 6.2pt)[d188864b5b69a3b51e3a66856a1db421941ac757a33304dbcfedbd6f1516116f]],
  [d064 replay `64..320`], [#text(font: "Menlo", size: 6.2pt)[2b11175e2cd775cac030fa5e38ef11f67b536a42a69a56bd01c646500f0da8a6]],
  [Certified segment `320..16,704`], [#text(font: "Menlo", size: 6.2pt)[514e8b0d6236ab883eb11f4cf24d8aeb66ce2a3ef6efb41652fcdd9a4fa08961]],
  [Generation receipt `320..16,704`], [#text(font: "Menlo", size: 6.2pt)[99f09bdad81ff06641063e1cdbc4d3b106622c97d9685b30c218bf873889ea0f]],
  [d064 replay `320..16,704`], [#text(font: "Menlo", size: 6.2pt)[93a36eeccdf6b821d4ffd8bebd552975ea2e0daa826672b2eb0dc88f95751dd8]],
  [Expanded certified ledger], [#text(font: "Menlo", size: 6.2pt)[822f6b0908432279351791e717dbc5caaf9600eb58bddc611834c63205dff5c2]],
  [Expanded uniqueness receipt], [#text(font: "Menlo", size: 6.2pt)[f014b8202ad435367824a5960387676bb3d568015e5e864f58cdde8c2bc069b1]],
  [Expanded sorted index], [#text(font: "Menlo", size: 6.2pt)[055fb9e007b0b7a261f4850b78450d1d8f306ec56e13c316cebf20514d1793ad]],
  [Expanded-prefix manifest], [#text(font: "Menlo", size: 6.2pt)[e26e618df5f47feabdf50534eccbe6fec7e1ba2b20581fedffec3f8a5b22a1b4]],
  [False-unconditional control], [#text(font: "Menlo", size: 6.2pt)[a8d615a2f5c8a2099d788820780c33531e44b05e8a5a81e8c357d57cf4227417]],
)

#v(10pt)
#text(size: 7.5pt, fill: slate)[The Markdown report in the same directory is the source-linked narrative. NAS paths and source identities are recorded there in full.]
