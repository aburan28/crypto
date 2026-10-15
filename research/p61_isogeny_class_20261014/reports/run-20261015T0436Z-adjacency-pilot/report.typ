#let navy = rgb("#17324d")
#let blue = rgb("#2367a8")
#let green = rgb("#15803d")
#let amber = rgb("#b45309")
#let slate = rgb("#475569")
#let pale = rgb("#f3f6fa")
#let rulecolor = rgb("#cbd5e1")

#set document(title: "p61 selected-degree adjacency", author: "Cryptanalysis research pipeline")
#set page(
  paper: "us-letter",
  margin: (x: 0.68in, top: 0.66in, bottom: 0.64in),
  header: context {
    if counter(page).get().first() > 1 {
      set text(size: 7.5pt, fill: slate)
      grid(columns: (1fr, auto), [p61 selected-degree adjacency], [Evidence report - 2026-10-15])
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
#set par(leading: 0.64em, justify: true)
#set heading(numbering: "1.")
#show heading.where(level: 1): it => block(above: 13pt, below: 7pt)[#set text(size: 15pt, weight: "bold"); #it]
#show heading.where(level: 2): it => block(above: 10pt, below: 5pt)[#set text(size: 11.5pt, weight: "bold"); #it]

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
  #v(0.22in)
  #text(size: 24pt, weight: "bold", fill: navy)[p61 selected-degree adjacency]
  #v(7pt)
  #text(size: 13.5pt, fill: blue)[Replayed production calibration and 6 TB corpus contract]
  #v(18pt)
  #rect(width: 84%, height: 1.4pt, fill: blue)
  #v(16pt)
  #grid(
    columns: (auto, 1fr),
    column-gutter: 14pt,
    row-gutter: 7pt,
    align: (right, left),
    [#text(weight: "bold", fill: slate)[Field]], [`p = 2^61 - 1 = 2305843009213693951`],
    [#text(weight: "bold", fill: slate)[Selected class]], [`p61-i0025`, trace `-704859453`],
    [#text(weight: "bold", fill: slate)[Class order]], [`358094304` vertices],
    [#text(weight: "bold", fill: slate)[Frozen degrees]], [`5, 13, 37, 59`],
    [#text(weight: "bold", fill: slate)[Run]], [`p61-adjacency-pilot-20261015T0436Z`],
  )
  #v(20pt)
  #tag("PRODUCER PASS", green)
  #h(7pt)
  #tag("INDEPENDENT REPLAY PASS", green)
  #h(7pt)
  #tag("FULL GRAPH ACTIVE / PARTIAL", amber)
]

#v(23pt)
#callout(
  "Decision",
  [A 4,096-vertex production window is admitted with 28,672 directed edges, 471,040 kernel coefficients, zero generation failures, and exact kernel-first replay. The calibrated representation projects to 3.632 TB logical JSONL and 0.748 TB gzip for the complete selected graph. Coverage is still 4,096 of 358,094,304 source vertices.],
  color: navy,
)

#v(10pt)
#align(center)[#text(size: 8pt, fill: slate)[Producer `f9b48fcf`; verifier `92a62be6`; different source commits, binaries, and reconstruction paths.]]

#pagebreak()

= Evidence state

#figure(
  image("evidence-pipeline.svg", width: 100%),
  caption: [Downstream graph and ML products remain gated on complete upstream evidence.],
)

#table(
  columns: (1.35fr, 0.82fr, 2.1fr),
  align: (left, left, left),
  table.header(
    table.cell(fill: navy)[#text(fill: white, weight: "bold")[Requirement]],
    table.cell(fill: navy)[#text(fill: white, weight: "bold")[State]],
    table.cell(fill: navy)[#text(fill: white, weight: "bold")[Boundary]],
  ),
  [Class size], [#text(fill: green, weight: "bold")[Certified]], [358,094,304 vertices; receipt bound into every shard],
  [Four-degree edge path], [#text(fill: green, weight: "bold")[Pass]], [1/2/2/2 rational edges per vertex on every admitted window],
  [Independent replay], [#text(fill: green, weight: "bold")[Pass]], [Kernel-first verifier at a distinct commit and binary digest],
  [Vertex census], [#text(fill: amber, weight: "bold")[Partial]], [2,659,968 records replayed through `x=147776`],
  [Global uniqueness], [#text(fill: amber, weight: "bold")[Partial]], [595,584 records externally deduplicated through `x=33088`],
  [Production adjacency], [#text(fill: amber, weight: "bold")[Partial]], [4,096 / 358,094,304 source vertices],
  [Dense join + reciprocity], [Pending], [Requires the complete globally unique vertex index],
  [ML-ready corpus], [Pending], [Labels, masks, split map, null controls, and joined targets remain],
)

= Frozen graph arithmetic

#table(
  columns: (0.7fr, 1.2fr, 1.1fr, 1.25fr),
  align: (right, left, right, right),
  table.header(
    table.cell(fill: pale)[#text(weight: "bold")[Degree]],
    table.cell(fill: pale)[#text(weight: "bold")[Behavior]],
    table.cell(fill: pale)[#text(weight: "bold")[Neighbors]],
    table.cell(fill: pale)[#text(weight: "bold")[Full edges]],
  ),
  [5], [ramified], [1], [358,094,304],
  [13], [split], [2], [716,188,608],
  [37], [split], [2], [716,188,608],
  [59], [split], [2], [716,188,608],
  [#text(weight: "bold")[Total]], [], [#text(weight: "bold")[7]], [#text(weight: "bold")[2,506,660,128]],
)

The complete selected graph also contains exactly 41,180,844,960 kernel coefficients. Other degrees, trace classes, and the field-wide curve population are separate scopes.

#pagebreak()

= Replayed calibration

#table(
  columns: (1.1fr, 0.72fr, 0.8fr, 0.86fr, 0.86fr, 1fr, 0.9fr),
  align: (left, right, right, right, right, right, right),
  table.header(
    table.cell(fill: navy)[#text(fill: white, weight: "bold")[Window]],
    table.cell(fill: navy)[#text(fill: white, weight: "bold")[Vertices]],
    table.cell(fill: navy)[#text(fill: white, weight: "bold")[Edges]],
    table.cell(fill: navy)[#text(fill: white, weight: "bold")[Generate]],
    table.cell(fill: navy)[#text(fill: white, weight: "bold")[Replay]],
    table.cell(fill: navy)[#text(fill: white, weight: "bold")[Logical]],
    table.cell(fill: navy)[#text(fill: white, weight: "bold")[Gzip]],
  ),
  [Path smoke], [1], [7], [18.703 s], [21.446 s], [12,227], [3,593],
  [Admission], [64], [448], [24.347 s], [23.712 s], [649,507], [135,394],
  [Production], [4,096], [28,672], [191.551 s], [296.931 s], [41,543,871], [8,552,556],
)

#v(8pt)
#figure(
  image("adjacency-calibration.svg", width: 100%),
  caption: [Single startup-inclusive observations on a shared host. Every window passed exact independent replay.],
)

The production window recorded 21.383 generation vertices/s, 13.794 replay vertices/s, and 8.385 vertices/s for the sequential admitted pipeline. A linear single-worker diagnostic is 193.8 generation days plus 300.5 replay days. It is not a schedule commitment: the eight-thread vertex census was active and no isolation or repeated timing design was used.

= Storage contract

#table(
  columns: (1.8fr, 1.15fr, 0.72fr),
  align: (left, right, right),
  table.header(
    table.cell(fill: pale)[#text(weight: "bold")[Projected product]],
    table.cell(fill: pale)[#text(weight: "bold")[Bytes]],
    table.cell(fill: pale)[#text(weight: "bold")[TB]],
  ),
  [Complete adjacency JSONL], [3,631,988,176,565], [3.632],
  [Complete adjacency gzip], [747,710,348,692], [0.748],
  [Logical allowance after adjacency], [2,368,011,823,435], [2.368],
)

The remaining allowance is reserved for useful dense joins, reciprocity audits, frozen attack contexts, observation masks, split maps, and ML-ready tables. Logical and compressed sizes remain separate. The supervisor enforces a 1.5 TB free-space floor and a 6.0 TB logical cap.

#pagebreak()

= Completion boundary

The vertex stream has per-segment replay through `x=147776`, but global uniqueness is sealed only through `x=33088`. The adjacency stream therefore records canonical target j-invariants without assigning dense target indices. The early `x=0..16704` prefix also requires a normalized source catalog or bounded backfill because it predates the campaign-admission schema.

Terminal completion requires:

+ every one of 358,094,304 sources exactly once;
+ exactly 2,506,660,128 directed edges and the fixed per-degree totals;
+ a globally unique dense vertex index with every target resolved;
+ reciprocal endpoint and degree consistency;
+ immutable manifests for every shard, retry, join partition, and exception; and
+ frozen hardness labels, explicit censor masks, field/class split maps, abstention, and null controls before M7 successor training.

= Preserved failures and limitations

- An exact-producer reconstruction under `/private/tmp` exhausted local scratch during archive extraction. The corrected detached worktree on SSD990-2 built in 602 seconds and reproduced the pilot binary digest exactly.
- The first supervisor unit-test invocation failed during dynamic module registration. The corrected final suite passes five tests, including the bounded-follow guard.
- The first persistent Terminal launch truncated the long inline command and exited before creating a log or child. A committed short wrapper corrected the interface failure; the supervisor, `caffeinate` guard, and second production-window child were alive at the recorded health observation.
- SMB per-file sync is unsupported and remains recorded. Host sync and post-write digest replay precede admission.
- Timings are descriptive shared-host observations, not a comparative speed result.

= Evidence index

#table(
  columns: (1.5fr, 2.5fr),
  align: (left, left),
  table.header(
    table.cell(fill: pale)[#text(weight: "bold")[Artifact]],
    table.cell(fill: pale)[#text(weight: "bold")[SHA-256]],
  ),
  [Producer binary], [`a69def0c21dc455077abb8bef86d828c13dfcbea8f45382f71614fdf260d7c78`],
  [Verifier binary], [`85b1265b187ecd2bbd2ff7fd0f7c9608d3ccd2e4a7234f73ee1b2863be599208`],
  [64-vertex admission], [`26b4358961f6cf3b58767b454a5272b101ea520a14ebb95ad409ded7bc37685f`],
  [4,096-vertex shard], [`ce0fd43d317b56e386462b7cb8ebd5605bba1195efaa57ba1e36819c5917bbc1`],
  [4,096 generation receipt], [`81fc61cf45c2718acd74236ebff5201f1e751c570c57afc86cdd1491a0526cbf`],
  [4,096 replay receipt], [`5755b14028bcf7bf51e19870ae3fa8e6422661e9836a6fc8780ead7e8878870e`],
  [4,096 admission], [`be33159cc6fb7bc3460c6626c827d2e01ea2cf46ada017fb8acafc19b21e26e8`],
)

#v(10pt)
#callout(
  "Durable root",
  [`/Volumes/data/ec-hardness-atlas-20261012/p61-selected-adjacency-v1/`],
  color: blue,
)
