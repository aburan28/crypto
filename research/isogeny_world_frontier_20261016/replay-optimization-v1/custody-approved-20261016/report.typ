#set page(paper: "us-letter", margin: (x: 0.64in, y: 0.55in))
#set text(size: 9.1pt)
#set par(leading: 0.35em)
#align(right)[TASK-20261016-de577e · 2026-10-16 · engineering checkpoint]
= Approved native custody
All four frozen comparison executables and the build/control inputs are preserved remotely and independently checked after full readback. The largest-prime-degree goal remains active; this checkpoint establishes no new degree or performance result.
#image("custody.svg", width: 100%)
#table(columns: (2fr, 3fr), inset: 5pt, stroke: 0.4pt + rgb("#aebfcd"),
 [*Requirement*], [*Recorded result*],
 [Initial logical names / distinct objects], [9 / 8; two dataset names share one content object],
 [Compiled executables], [4 uploaded; installed = compiled before and after controls],
 [Native remote size/checksum], [PASS for all 9 names],
 [Independent complete readback], [SHA-256 and byte equality PASS for all 8 distinct objects],
 [Build/control evidence], [61-file archive plus manifest and full accepted-P192/13 preparation datasets],
 [Scientific trials / new scores], [0 / 0],
 [Independent source review], [Not started; source gate remains incomplete],
)
*Destination:* #raw("s3://crypto-autoresearcher/store"). Logical prefix:
#raw("evidence/isogeny-fixed-modulus-20261016/015fa7bb504e/")

*Provenance.* Frozen source: #raw("015fa7bb504ef63953e0e04e464bb266f4f207bf"). Exact binary identities, physical S3 URIs, byte sizes and readback paths are in verification.json; inputs.json binds the original files. The successful native transfer exited 0 at 18:57:15 UTC. Explicit user approval preceded it. The earlier denied attempt executed zero uploads.

*Audit.* Native push receipts, complete fetch/verify stdout and stderr, and the store index are retained. Python hashlib and direct byte comparisons independently checked downloads. The native transport executable is hash-bound; its build-source identity was not independently established. Python performed metadata/hash/packaging only. No endpoint override or replacement writes were used.

The initial metadata assembler incorrectly required stderr to be empty; it stopped on ordinary native fetch progress after byte/hash checks. Its source and the corrected successful checker are both retained. Native transfer itself was not retried. Separate administrative checks encountered missing PyYAML and system-temp exhaustion; the existing project virtualenv and metadata commands that avoid here-document temporaries resolved them without deleting evidence; scratch remains on SSD990-2. Pull/verify plain-text files were named .stdout before sealing. Read README.md and metadata-assembly-attempts.md for exact reproduction and this failed assembly attempt.
#pagebreak()
= Windows fixture and current-source gate
Windows CI run 38933315615 / job 117082574923 failed before mathematical replay at the accepted P224/1471 fixture-hash assertion. Actual Git checkout with autocrlf=true and eol=crlf reproduced the observed CRLF digest. Marking exactly three fixture paths as -text made fresh checkouts match the original Git blobs byte for byte.

Expected LF digest: #linebreak()
#text(size: 7.6pt)[#raw("ecf1109dbff2a72000fd84d9939127db68e2287e55e4ec10d143229874c230c8")]

Observed and reproduced CRLF digest: #linebreak()
#text(size: 7.6pt)[#raw("854333812fd0dffb0bb3479d545d36f7805ee39c85314da57d6db59607f7fb5c")]

The expected hash, fixture bytes, mathematical code, assertions and admission requirements are unchanged. A Windows-host retest is pending.
#table(columns: (1.3fr, 2.8fr), inset: 5pt, stroke: 0.4pt + rgb("#aebfcd"),
 [*Frozen native arm*], [*Existing cache/adversarial integration test*],
 [Reference], [1 passed / 0 failed / 19 filtered],
 [Candidate], [1 passed / 0 failed / 19 filtered],
 [Full verification output], [Byte-identical JSON; 0 exclusions; all 20 public transports],
 [Raw evidence], [106 files in checked archive; hashes in native-manifest.json],
)
These are existing native executable tests on a shared Apple host. Cache activation and both rehashed negatives are exercised. The separate control archive also has verified remote custody; its store index and full readback proof are in ../windows-fixture-20261016/custody/. All four compiled/installed executable hashes remain unchanged. Operational OS receipts are retained; no timing ratio or isolated performance result is inferred.

*Current source.* The checkpoint checkout merged main at 3ce9694072773badcedc815ddbb953339888e433. source-delta.json binds differences from frozen 015 source: separate CM/CRT API/CLI/test/doc additions, with dedicated replay arithmetic and previous API bodies unchanged. The root library also received unrelated hwcost/inherited-F4/Koblitz changes, whose fresh complete serial release suite finished 4161 passed / 1 failed / 117 ignored, zero filtered, exit 101. The sole SQISign statistical assertion observed 27/4000 = 0.00675 below its strict lower bound 1/(37*4). Its module and SHA3 dependency are unchanged against frozen source and the main parent; cause remains unestablished. One predeclared whole-module retake on the same executable passed 14/0/0, 4265 filtered. It cannot replace the failed full suite. root-gate.json and its remotely verified raw archive retain both runs, source comparisons and OS receipts. Administrative task-release tests passed 24/24 at separate HA source69928.

*Preserved history.* Prior root default-parallel run: 4146 passed / 2 failed / 117 ignored. Complete serial retake: 4148 / 0 / 117; focused JV group: 14 / 0 / 1. Variation's cause and parallel stability are not established. IDC1h41e9973f839c85bf and two checked Lean lemmas remain bounded; formal-series/reversal bridges remain open.

*Next.* Coordinator snapshot c447e0, then fresh independent source/resource review 06cc13 with prospective source/raw-first independence. Producer conclusions must not seed that review. Source gate, later protocol approval and board promotions remain Coordinator decisions. No 50341 replay, larger construction, compact target or measurement campaign ran here. All 54 standard prime-field boards and canonical curve/cost graphs remain unchanged; no new graphable mathematical result was established. Full evidence: README.md, verification.json, source-delta.json, root-gate.json and ../windows-fixture-20261016/. The named problem records the unresolved root gate; its row is included in the shared problem register after that scope became available. All 11 logical names / 10 distinct custody objects passed full readback checks. Producer delivery does not establish source admission or a green PR.
