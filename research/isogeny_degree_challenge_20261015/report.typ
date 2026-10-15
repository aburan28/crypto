#set page(paper: "us-letter", margin: 0.7in, numbering: "1")
#set text(font: "DejaVu Sans", size: 10pt)
#set par(justify: false)
#let run = sys.inputs.at("run")
#let data = json(run + "/audit.json")
#let freeze = json(run + "/freeze.json")
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

#figure(image(run + "/admission.svg", width: 100%),
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

= Reproduction and evidence
Use the dedicated ReleaseMe #raw("isogeny") executable. The runtime embeds its
catalogue, rules, constructor and replay engines and requires no Cargo checkout.
Publication of this new command remains subject to the main release workflow.

#raw("isogeny challenge policy\nisogeny challenge plan --curve all --from 10454 --to 65537\nisogeny challenge run --curve p192 --ell 13 --out fresh-control\nisogeny challenge judge --input submission --out fresh-judge\nisogeny challenge leaderboard --input submissions.json --out fresh-board", block: true)

Source-linked evidence lives beside this report: audit.json, freeze.json,
policy.json, plan.json, raw native control packets, the fresh historical-anchor
replay and validation logs. The Markdown report identifies exact archive paths
and test outcomes. audit.rs reproduces coverage and rejection counts.

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
