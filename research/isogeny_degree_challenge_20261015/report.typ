#set page(paper: "us-letter", margin: 0.7in, numbering: "1")
#set text(font: "Libertinus Serif", size: 10pt)
#set par(justify: false)
#let run = sys.inputs.at("run")
#let data = json(run + "/audit.json")
#let freeze = json(run + "/freeze.json")
#let wrapped_id(value) = {
  for start in range(0, value.len(), step: 36) {
    text(size: 8pt, value.slice(start, calc.min(start + 36, value.len())))
    if start + 36 < value.len() { linebreak() }
  }
}
#align(center)[#text(size: 22pt, weight: "bold")[Prime Degree Frontier]\
Native challenge and admission controls · 15 October 2026]

= Scope and frozen identity
The challenge ranks the largest independently verified odd prime degree
separately for every exact registered standard prime-field source curve. Expanded
normalized rational maps and monic kernel polynomials are required. Each
leaderboard build runs fresh independent replay; saved scores never enter ranking.

#table(columns: (1.4fr, 4fr), inset: 5pt,
  [Source revision], [#text(size: 8pt)[#freeze.source_commit]],
  [Rules/catalogue season], [#text(size: 8pt)[#data.season]],
  [Executable SHA256], [#text(size: 8pt)[#freeze.binary_sha256]],
  [Run host], [#freeze.host],
  [Run scope], [Correctness/admission controls on a shared host; no throughput comparison.],
)
#link("https://github.com/aburan28/crypto/blob/" + freeze.source_commit + "/tools/isogeny-cli/src/challenge.rs")[Frozen judge source] · #link("audit.json")[Native audit] · #link("raw.tar.gz")[Raw packet]

= Observed coverage and controls
#table(columns: (3fr, 1fr), inset: 5pt,
  [Registered standard prime models], [#data.registered_prime_models],
  [Construction adapter available], [#data.construction_models],
  [Independent replay available], [#data.replay_models],
  [Replay adapter unavailable, retained in roster], [#data.unavailable_replay_models],
  [Candidate curve-degree pairs above 10453 through 65537], [#data.shortlisted_candidates],
  [Altered-artifact controls rejected], [#data.rejected_mutations.len()],
  [Verified source rows in control leaderboard], [#data.verified_source_rows],
  [Short-map scheduling controls accepted / executed], [#data.resource_stress_accepted / #data.resource_stress_count],
  [Short-map scheduling controls censored], [#data.resource_stress_censored],
)
Counts are exact native auditor outputs, not sampled estimates. The low-order
Frobenius screen produces candidates only. Missing adapters remain explicit.
Two known-answer fixtures score P192/13 and P521/7 on different source rows.
These controls contribute zero new discoveries. A forged saved score is rejected.

#for control in data.known_answer_controls [
  #text(size: 8pt)[Exact source: #control.source_icv1; verified degree #control.prime_degree.]\
]

#figure(image(run + "/admission.svg", width: 78%),
  caption: [Native admission-control flow. Roster and control counts come from audit.json.])

#pagebreak()
= Rejection evidence and resource policy
#table(columns: (1.2fr, 1fr), inset: 5pt,
  [Control], [Observed result],
  ..data.rejected_mutations.map(r => (r.control, r.status)).flatten(),
)
Kernel, numerator and source mutations are tested after recomputing their
receipt hashes. Other controls alter the degree, substitute a compact-only
claim, traverse a curve path, or omit a hash. A prewritten PASS and an inflated
stored score do not change the freshly verified degree. Output-directory
reuse is refused and the first judgement's bytes remain unchanged.

Replay has a 3600-second wall budget and an 8 GiB sampled RSS cap per submission.
A lower replay timeout is permitted. The bundled constructor has a separate
1800-second/8 GiB limit. The raw map limit is 256 MiB, summary 1 MiB and stderr
8 MiB. Resource failures are censored with no score; invalid maps are rejected.
Snapshot copying precedes the bounded child process. Operational receipts are
retained; no isolated timing or peak-memory guarantee is inferred.

The bundled constructor accepts degrees through 1000000. The judge accepts
u64 prime degrees only when the full explicit artifact fits the byte and
verification limits. These implementation limits do not prove nonexistence.
External constructor costs are not independently attested or used as a tiebreak.

#if data.starting_curve_rows.len() > 0 [
  #pagebreak()
  = Fresh starting leaderboard
  These historical maps were replayed again under this run's pinned season.
  Each degree scores only its own exact source. The short-budget P192/10453
  control returned #data.short_budget_control.status with no score; the
  full-budget leaderboard used a separate fresh replay.
  #table(columns: (3fr, 0.7fr, 0.9fr, 0.5fr), inset: 5pt,
    [Exact source ICV1], [Prime degree], [Coverage], [Maps],
    ..data.starting_curve_rows.map(r => (wrapped_id(r.source_icv1), r.prime_degree, r.coverage, str(r.map_count))).flatten(),
  )
  Every admitted map passed the independent kernel, codomain and exact-map
  checks plus twenty public scalar transports. Re-admitting an existing map
  contributes zero additional discoveries. The other roster rows stay unscored.
  #link("raw.tar.gz")[Raw packet and fresh per-entry certificates]
  #parbreak()
  #for row in data.starting_curve_rows [
    #text(size: 8pt)[Degree #row.prime_degree: #row.source_icv1]
    #linebreak()
    #for target in row.targets [
      #text(size: 8pt)[→ #target]
      #linebreak()
    ]
  ]
]

= Reproduction and evidence
Use the dedicated ReleaseMe #raw("isogeny") executable. The runtime embeds its
catalogue, rules, constructor and replay engines and requires no Cargo checkout.
Publication of this new command remains subject to the main release workflow.

#raw("isogeny challenge policy\nisogeny challenge plan --curve all --from 10454 --to 65537\nisogeny challenge run --curve p192 --ell 13 --out fresh-control\nisogeny challenge judge --input submission --out fresh-judge\nisogeny challenge leaderboard --input submissions.json --out fresh-board", block: true)

Source-linked evidence lives beside this report: audit.json, freeze.json,
compressed policy/plan JSON and the complete raw archive. That archive contains
the native control packets, historical-anchor replay and attempt receipts.
The Markdown report identifies exact archive paths and test outcomes.
audit.rs reproduces coverage and rejection counts.

= Limits and related graphics
This is a challenge for pinned registered models, not a worldwide-record claim
or a finite mathematical maximum. Full explicit maps are its sole accepted
representation. Compact or custom-curve tracks require a separate contract.
The IC/rho scoreboard, IC leaderboard, progress timeline and canonical curve
edges remain unchanged: this packet establishes admission behavior and retains
known-answer maps without a new IC ratio, new registry model or new curve edge.
The earlier CPU optimization and additional-map campaigns remain separate.

The first full-suite attempt censored a valid short replay with zero RSS samples.
A later installed run on that source passed, confirming a scheduling-dependent
startup race. The fixed source makes the replay child wait until its supervisor
records positive RSS. The follow-up packet adds twenty fresh judgements of the
same public P192/13 map, each requiring a real sample. Earlier failed build/test
logs and the pre-fix packet remain archived; they are not rewritten as passes.
The intermediate stress packet accepted 19 of 20 and censored one during OS
exit cleanup. The final correction uses a completion marker after all algebra
and failure handling. Any missing-monitor refusal during active work stays
censored; a later passing certificate never erases it. Zero executed stress
controls in the table means that the earlier run did not request this phase.
