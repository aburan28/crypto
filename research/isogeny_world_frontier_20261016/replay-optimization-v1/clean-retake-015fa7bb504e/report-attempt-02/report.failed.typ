#set page(paper: "us-letter", margin: (x: 0.63in, y: 0.56in))
#set text(size: 9.1pt)
#set par(leading: 0.38em)
#align(right)[#text(size: 8pt)[TASK-20261016-de577e · clean source015fa7bb504e · 2026-10-16]]
= Clean native source/control checkpoint
The largest-degree goal remains active. This checkpoint preserves exact default/reference behavior and validates the default-off cached kernel-reciprocal implementation. It establishes no new admitted degree or performance result.

#image("../architecture.svg", width: 100%)

#table(columns: (2fr, 2.7fr), inset: 5pt, stroke: 0.4pt + rgb("#aebfcd"),
 [*Evidence gate*], [*Observed outcome*],
 [Clean build and installed byte identity], [Both policies dirty=false; compiled bytes equal installed bytes],
 [Complete native suite], [128 passed,0 failed per arm; one parent-invoked child helper ignored],
 [Complete mathematical output], [Six byte-identical JSON pairs; zero exclusions; all20 transports per map],
 [Rehashed/adversarial negatives], [Nine trace pairs plus two actual-cache negatives agree],
 [Lost/zero-RSS controls], [16 retained fail-closed injection receipts per arm; two worker shapes],
 [Deadline controls], [Two retained timeout receipts per arm],
 [Exact-source Lean], [Two checked lemmas; standard propext/Quot.sound only],
 [Remote custody / independent review], [LOCAL ONLY / NOT STARTED; source gate incomplete],
)

*Exact source:* 015fa7bb504ef63953e0e04e464bb266f4f207bf.
Both builds use rustc1.93.1, the same lockfile, target and release profile; only the one explicit feature differs. The source and version JSON identify the policy separately from mathematical certificates.

*Executable SHA256.* Reference:
#text(size: 7.8pt)[#raw("2ff7d57174cda2d0e90088eb6ce217cbd1826e7697ed12ec7770d3d8f395d379")]
Candidate:
#text(size: 7.8pt)[#raw("305bdff7db7bc25a2ccdce3f44e71c3d47269e76afa0ecd1703803365ab9b207")]

Native manifest and all control hashes: native-manifest.json. Native archive:
1,261 regular files,2,049,977 bytes; every archived hash replayed. Executables and target directories are excluded; the compiled Lean object is included. The archive and exact receipt are beside this report.
#pagebreak()
= Root failures, retakes and limits
The clean default-parallel root library run failed4146/2, with117 ignored. Both failures were unchanged JV solver admission-count tests (21 and22 admitted cases versus thresholds25 and30). The complete JV group then passed14/0/1 ignored serially. The entire root library then passed4148/0/117 ignored with one test thread, without skipping a test or changing any threshold or source.

The same root test executable was used by the earlier default-parallel pass, the later default-parallel failure and the full serial pass. Its SHA256 is
#text(size: 7.7pt)[#raw("e6b0366e26794fe703741ba309cd3c7cb168150c827fa2a23ba068038aec7399")]
Both the failed result and successful retakes are retained in root-gates/. This is observed variation; its cause is not established. Parallel stability is not asserted.

*Accounting.* The complete serial library body reported1047.40s. Its OS receipt recorded1047.87s real,799.80s user and101.98s system. The intervals overlap. Native builds/tests retain their own OS receipts. These are engineering costs on a shared Apple M4 Pro host, not isolated performance observations. Earlier development CPU costs remain incomplete.

*Input-preparation control.* Native generation used the accepted P192/13 map only, produced eight dataset cases, and returned PASS with zero unavailable RSS samples for both compiled variants. The two datasets are byte-identical:
#text(size: 7.7pt)[#raw("2b5bb3f2cdd04ddb3ab068e99c7a81e22e3801bcfd0471759ca993ec9e5fd871")]
Generation computes reference expectations in both variants. These bounded controls do not measure cached speed. The large datasets remain locally retained and hash-bound in generation-control.json.

*Proof scope.* IDC1h41e9973f839c85bf is a finite-field evaluation supplement. Lean checks the scalar integer identity and an integer coefficient-convolution precision lemma. Full formal-series uniqueness and quotient-reversal formalization remain open; no fresh lean4checker was available.

*Failures preserved.* The initial native alias failure, initial Lean lemma-name failure, later metadata path error, root parallel failure and all corrected controls remain accessible. The first native archive was sealed before the last full-root retake; that final retake is a separate immutable addition.

*Next gate.* Remote custody awaits authorization, and independent source/resource review has not started. TASK-de577e is not complete. No measure/analyze campaign,50341 extraction/replay, larger construction, compact target, default activation or board promotion ran. Source-only PR1991 remains the implementation PR. All 54 standard prime-field boards and the separate reviewed M44497 compact result are unchanged.

