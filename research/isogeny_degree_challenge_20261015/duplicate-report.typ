#set page(paper: "us-letter", margin: 0.7in, numbering: "1")
#set text(font: "Libertinus Serif", size: 10pt)
#let run = sys.inputs.at("run")
#let d = json(run + "/diagnostic.json")
#let entry = json(run + "/observed-entry.json")
#let cert = json(run + "/observed-certificate.json")
#align(center)[#text(size: 21pt, weight: "bold")[Prime Degree Frontier]\
Adversarial completeness regression]
= Observed outcome
The frozen ff71f5772 executable accepted a modified public P192/13 artifact
whose two map slots contained the same valid map. Fresh verification proved
the degree-13 map in each slot, but the COMPLETE coverage label was unsupported.
The native regression expected rejection and failed with observed exit 0 rather
than 1. This is a coverage-metadata defect; it does not invalidate the degree.
#table(columns: (1fr, 3fr), inset: 5pt,
 [Executed binary source], [#text(size: 8pt)[#d.binary_source]],
 [Test source], [#d.test_source],
 [Run time basis], [#d.observation_time_basis: #d.observed_run],
 [Fresh judge result], [#entry.status; degree #entry.prime_degree; coverage #entry.degree_coverage],
 [Declared maps / distinct kernels], [2 / 1],
)
= Input and controls
The native test copied isogenies[0] into isogenies[1] of the real public
P192/13 constructor output and recomputed its construction receipt SHA256.
The standard source, map coefficients, codomain and subgroup transports stayed
mathematically valid. The two displayed records therefore described one kernel
rather than two distinct Frobenius eigenlines.
#for record in cert.fresh_certificate.records [
 #text(size: 8pt)[Source #record.source_icv1 → #record.target_icv1]
 #linebreak()
]
#figure(image("duplicate-flow.svg", width: 100%), caption: [Observed false completeness claim. Both slots contained the same valid map.])
= Correction and acceptance
The corrected judge compares kernels modulo the source characteristic before
accepting two-line coverage, and requires two distinct Frobenius roots for a
COMPLETE claim. A partial valid map remains eligible; repeated-root totals are
left unestablished. The final corrected-source packet must pass the original
controls plus this ninth negative case, twenty scheduling controls and fresh
starting-map replay. No additional discovery or speed result is credited.

The exact altered bytes, updated receipt, fresh certificate, accepted entry,
stdout/stderr and failed assertion are in #link("raw.tar.gz")[the raw packet].
#link("observed-entry.json")[Observed entry] · #link("observed-certificate.json")[Fresh certificate]
The code and diagnostic use public known-answer fixtures on a shared M4 Pro;
there is no isolated timing comparison. The root IC/rho and curve-edge charts
remain unchanged because this is an admission-correctness correction.
