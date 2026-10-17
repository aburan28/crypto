# WS6 public literature and incident corpus protocol

> **Status: draft, non-authoritative, not frozen.** Protocol revision
> `ws6-corpus-protocol-draft-2026-10-13`; prepared 2026-10-13. No search has
> been executed under this revision. It establishes neither corpus
> completeness nor the truth of any seeded claim. `WS6-G1` and `WS6-G2`
> remain open.

This protocol preregisters how Workstream 6 will discover, acquire, classify,
and review public evidence about elliptic-curve discrete-logarithm methods and
other elliptic-curve private-key recovery routes. It covers public literature,
standards, source code, advisories, and public incident records. It does not
authorize live-target testing, acquisition of private incident data, or use of
recovered secrets.

The companion files are:

- [`README.md`](README.md), the authority boundary, artifact index, and local
  syntax/integrity checks;
- [`source-claim-manifest.schema.json`](source-claim-manifest.schema.json), a
  strict structural contract whose draft authority flags fail closed;
- [`seed-manifest.json`](seed-manifest.json), five pointer-only examples drawn
  from citations already present at repository revision
  `1e447956d3616092a223264f1c6c31eff8ac9979`; and
- [`fixtures/`](fixtures/), an empty valid draft and four schema-negative
  controls for authority, completeness, freeze, and review gates.

The seed is a schema fixture and discovery queue, not an executed search, an
endorsed bibliography, an incident catalogue, or a comparator matrix.

## Requirement-to-evidence status

| Requirement | Artifact in this draft | Status | What remains |
|:--|:--|:--|:--|
| Preregister databases, queries, dates, languages, screening, deduplication, and citation passes | This protocol | **implemented but unverified** | freeze an execution revision, review it, and retain run receipts |
| Bind sources and claims with exact provenance and fail-closed status | Draft JSON Schema and fixtures | **implemented but unverified** | independent schema and semantic review; native verifier |
| Seed only from checked-in primary-source pointers | Seed manifest, pinned to three Git blobs | **verified for the five seed pointers** | retrieve and authenticate the actual sources; the claims remain `unknown` |
| Build the systematic source and incident corpus | none | **not attempted** | execute all frozen searches and citation passes |
| Establish comprehensive coverage | none | **not established** | freeze denominators and comparator versions, resolve gaps, and obtain second-person review |

## 1. Questions and classification boundary

The corpus answers evidence questions, not whether a curve or implementation is
safe:

1. What public primary sources state a falsifiable claim relevant to an audit
   route in lanes A--E of the master plan?
2. What exact mathematical, curve, subgroup, protocol, implementation, version,
   access, and cost hypotheses delimit that claim?
3. Which later source corrects, refutes, supersedes, reproduces, or narrows it?
4. For a public incident, what is observed by the affected party or maintainer,
   what is inferred by a third party, and what remains unattributed?
5. Which enumerated query results and comparator entries have not yet been
   acquired or screened?

Every record must distinguish these route lanes:

- `A`: an algorithm for the original ECDLP;
- `B`: a mathematical transfer, descent, cover, pairing, isogeny, or lift;
- `C`: a parameter or structural precondition;
- `D`: a protocol or implementation bypass, including nonce, validation,
  side-channel, fault, RNG, firmware, and supply-path failures; and
- `E`: an evidence, identity, provenance, custody, build, or accounting issue.

Recovering a private key through lane D is not an ECDLP solution. Existence of
a lane-B map is not evidence that it preserves the target subgroup or gives a
cheaper destination solver. A failed bounded search is evidence only for its
recorded query, date, database snapshot, language, and screening rules.

## 2. Freeze before execution

The first executable search PR must copy this draft to a versioned protocol and
freeze the following before issuing any query:

- protocol identifier, repository commit, reviewer, and UTC freeze time;
- each service, endpoint or user interface, database release/snapshot where
  available, access tier, and query syntax version;
- the exact query string sent, field restrictions, filters, sort order,
  pagination rule, and result cap;
- inclusive coverage dates and a single UTC evidence cut-off;
- language set and translation procedure;
- inclusion, exclusion, deduplication, and reviewer-conflict rules;
- backward/forward citation-pass budgets and saturation condition;
- comparator names, versions or retrieval timestamps; and
- storage, lawful-access, content-hash, and redaction rules.

Any post-freeze change creates a new protocol revision. It must state the reason
and re-run the affected cells; it must not silently alter the denominator. A
pilot used to debug syntax is logged as `pilot` and cannot contribute included
records unless repeated under the frozen protocol.

## 3. Planned public sources

The table is a preregistered source registry, not proof that access exists or
that a search has run. The executable revision must pin the actual endpoint and
version for every used row and record unavailable rows as gaps.

| Class | Planned services or collections | Search role | May alone support a critical claim? |
|:--|:--|:--|:--|
| Cryptographic and mathematical literature | IACR ePrint, DBLP, Crossref, OpenAlex, IEEE Xplore, ACM Digital Library, SpringerLink, zbMATH/OpenAlex mathematics metadata | discovery, metadata, forward/backward chaining | only the retrieved primary paper, formal erratum, or author version; index metadata alone cannot |
| Standards and guidance | IETF Datatracker/RFC Editor, NIST CSRC, SECG/Certicom publications, relevant BSI and national standards portals | requirements, parameter provenance, revisions, errata | yes, for what that exact version normatively states |
| Vulnerabilities and public incidents | CVE records and CNA references, affected-vendor advisories, maintainer releases/commits, CERT/CC, CISA KEV, NVD | identifier discovery, affected-version bounds, public incident evidence | vendor/CNA/maintainer primary evidence may; aggregators alone cannot establish impact or attribution |
| Code and artifacts | official source repositories, signed releases, commit histories, issue trackers, reproducibility archives | implementation identity, fix provenance, executable evidence | only when exact revision, artifact identity, and claim locator are bound |
| Comparators | frozen SafeCurves pages, ECRYPT/ENISA guidance, major public ECC bibliographies, standards catalogues, and public advisory corpora | gap finding and denominator comparison | no; follow their citations to a primary source |

Discovery sources and comparators are retained with `source_role` set to
`secondary` or `comparator`. They cannot be silently relabelled primary.

## 4. Preregistered query families

Each service receives a syntax-correct expansion of every applicable family.
The exact expansion, including escaping and filters, is stored in the search-run
receipt. Quoted strings below are conceptual templates, not claims that all
services accept the same grammar.

| ID | Conceptual query template | Primary lanes |
|:--|:--|:--|
| `Q-A-GENERIC` | `("elliptic curve discrete logarithm" OR ECDLP) AND (algorithm OR attack OR complexity OR Pollard OR kangaroo OR "index calculus" OR "summation polynomial")` | A |
| `Q-B-TRANSFER` | `("elliptic curve" OR ECDLP) AND (MOV OR "Frey Rueck" OR "Frey Ruck" OR pairing OR anomalous OR p-adic OR "Weil descent" OR GHS OR isogeny OR twist OR cover OR hyperelliptic OR Jacobian OR restriction-of-scalars)` | B, C |
| `Q-C-PARAMETER` | `("elliptic curve" OR ECDSA OR ECDH) AND (embedding-degree OR trace-one OR anomalous OR supersingular OR cofactors OR twist-security OR endomorphism OR parameter-generation)` | C |
| `Q-D-NONCE-RNG` | `(ECDSA OR EdDSA OR "elliptic curve signature") AND (nonce-reuse OR nonce-bias OR "hidden number" OR RNG OR entropy OR deterministic-nonce OR lattice)` | D |
| `Q-D-VALIDATION` | `(ECDH OR ECDSA OR X25519 OR Ed25519 OR "elliptic curve") AND (invalid-curve OR twist OR subgroup OR cofactor OR validation OR parser OR "signature verification")` | D |
| `Q-D-PHYSICAL` | `(ECDSA OR ECDH OR "elliptic curve") AND (timing OR cache OR side-channel OR power OR EM OR electromagnetic OR fault OR firmware OR HSM OR TPM)` | D |
| `Q-D-INCIDENT` | `(ECDSA OR ECDH OR "elliptic curve") AND (CVE OR advisory OR vulnerability OR incident OR compromise OR "key recovery")` | D, E |
| `Q-E-EVIDENCE` | `("elliptic curve" OR ECDLP) AND (erratum OR retraction OR reproduction OR benchmark OR artifact OR provenance OR "supply chain")` | E |

For advisory systems that support structured fields, execute parallel searches
over product, weakness, keyword, and reference fields instead of relying only
on free text. For author/title indexes, include exact-title and cited-author
queries derived from included sources. Synonyms added after screening begins are
a protocol amendment and trigger the affected searches again.

## 5. Dates and languages

- The planned lower date bound is the earliest indexed record; no historical
  lower cut-off is imposed.
- The draft evidence cut-off is `2026-10-13T23:59:59Z`. The executable revision
  must confirm or replace it before freezing. Retrieval after that time may
  authenticate an in-window record but must not add later publications to that
  revision.
- Initial title/abstract screening languages are English, French, German,
  Spanish, Chinese, Japanese, and Russian. A record is not excluded solely for
  language. If competent review or a traceable translation is unavailable, use
  exclusion reason `deferred_language_review`, retain the pointer, and count it
  as an unresolved gap.
- Queries must use translated route terms where the service has meaningful
  non-English coverage. Machine translation may assist discovery but may not be
  quoted as the primary source or used alone to verify a critical claim.

No global completeness percentage may be reported from this language set.

## 6. Screening and eligibility

Two passes are retained separately: title/metadata screening and full-source
screening. Each decision records reviewer, UTC time, protocol revision, reason
code, and source identifier. A random 10% of metadata exclusions (minimum 25 if
available) and every critical-claim inclusion receive second-person review.
Disagreements remain open until reconciled with a recorded reason; raw decisions
are never overwritten.

Include a source when it is public and at least one of the following holds:

- it states or proves a mathematical route, bound, applicability condition, or
  counterexample relevant to lanes A--C;
- it specifies elliptic-curve parameters, validation, nonce, key-generation,
  or protocol behavior needed to interpret a route;
- it is primary code, a release, fix, advisory, CNA record, or affected-party
  incident statement relevant to lanes D--E;
- it is a formal erratum, retraction, correction, reproduction, or conflicting
  result for an included claim; or
- it is a frozen comparator used to audit search coverage.

Exclude from the claim-supporting primary corpus, while retaining the discovery
record and reason, when it is:

- unrelated after full-source review (`out_of_scope`);
- duplicate bytes or the same version under another URL (`duplicate_exact`);
- a superseded copy retained under its successor (`superseded_version`);
- commentary with no independently locatable primary support
  (`secondary_only`);
- inaccessible within the frozen lawful-access budget (`unavailable` or
  `paywalled_unacquired`); or
- non-public or obtained under terms that do not permit this corpus record
  (`nonpublic_or_restricted`).

Unavailable sources remain counted gaps. They are not evidence for or against
their alleged claim.

## 7. Source identity, acquisition, and custody

Resolve source identity in this order:

1. primary stable identifier and version (DOI, IACR ePrint revision, RFC plus
   update/errata state, standard plus edition, CVE plus CNA record, repository
   commit, or signed release);
2. lawful retrieved-byte SHA-256, byte count, media type, and retrieval time;
3. issuer, title, creators, publication date, and version metadata; then
4. URL only as a locator, never as the sole identity if a stronger identifier
   exists.

Exact duplicate bytes merge only at the acquisition layer. Preprint, conference,
journal, erratum, revised standard, advisory update, patched release, and
retrospective incident report remain separate source records linked by
`supersedes`, `corrects`, or `is_version_of`. A withdrawn or retracted source is
retained with that status.

For each lawful acquisition, record append-only custody events: request URI or
repository identity, retrieval method, actor or automation identity, UTC time,
byte count, media type, SHA-256, license/access observation, and any
normalization. Hash the original bytes. If only metadata or a checked-in pointer
is available, use `pointer_only`; leave content hash, byte count, and retrieval
time null. Never invent a hash for a web page that was not captured.

## 8. Claims and exact locators

Sources and claims are many-to-many. A claim is an atomic, falsifiable sentence
with explicit route lane, hypotheses, affected identities, limitations, status,
and evidence maturity. Each source--claim link supplies one exact locator:

- page/theorem/equation/section for papers;
- section/requirement/table/erratum for standards;
- file and commit-bound line or symbol for code;
- advisory field, CNA record field, release note, or incident-report section;
  or
- `bibliographic_pointer_only` for an unacquired discovery seed.

Use `metadata_only` when no content was retrieved. Such a link cannot verify a
substantive claim. Quotations preserve the exact language and page; paraphrases
are labelled. Claims record `supported`, `refuted`, `unknown`, or
`not_applicable` separately from maturity `E0`--`E7`. A critical claim requires
the primary bytes or exact immutable source revision, an exact locator, and
second-person verification before it can be marked supported or refuted.

Every manifest binds the WS5 route registry separately from the coarse A--E
lane tags. Each source and claim carries a `route_registry_join` object. An
`unresolved` join has no route IDs, reviewers, or review time and remains
non-authoritative. A `candidate_unreviewed` join may name route IDs but cannot
be treated as resolved. A `reviewed_non_authoritative` join requires named
reviewers and a review time, but still cannot promote a claim. Before any join
is consumed, the semantic verifier must bind a SHA-256 digest of the registry's
raw file bytes (`digest_scope="raw_file_bytes"`), resolve every route ID, check
that the route lane agrees with the record's lane set, and retain the mapping
review receipt. The schema requires a hash-verified binding for every candidate
or reviewed join and a reviewed join for a supported or refuted claim. The seed
intentionally stops at a path-only registry binding and unresolved joins; it
does not invent route IDs.

The maturity labels refer to the non-authoritative WS4 definitions in
[`../evidence/evidence-profile-registry.draft.json`](../evidence/evidence-profile-registry.draft.json).
Maturity and truth status are orthogonal. In the seed, `E0` means only that a
discovery record and its checked-in pointer have been identified at a pinned
repository revision; it does not authenticate the cited source or support the
queued substantive claim. A later validator must reject any claimed maturity
whose WS4 gates are not actually evidenced.

The draft schema is intentionally non-authoritative. JSON Schema checks shape,
types, enums, and local fail-closed fields; it cannot prove that IDs resolve,
registry or source hashes match, route lanes agree, links are acyclic, counts
reconcile, or reviewers are independent. A native semantic verifier with
invalid fixtures is required before `WS6-G2` can close.

## 9. Citation passes and deduplication

For every included primary paper or incident source:

1. perform one backward pass over its references;
2. perform one forward pass using every frozen citation index that supports it;
3. screen corrections, errata, retractions, replications, and cited software;
4. repeat a second backward/forward pass only for newly included critical
   sources; and
5. stop at the preregistered pass limit, logging unexamined edges.

Saturation means a complete pass added no new *eligible records within the
frozen services and budgets*. It is not universal literature saturation.

Deduplicate first by exact stable identifier and version, then by content hash,
then by reviewed normalized metadata. Title similarity alone proposes a merge;
it cannot perform one. Incident aliases, CVE records, vendor advisories, fixes,
and public incident narratives remain distinct and are linked, because they may
make different claims and change at different times.

## 10. Search-run receipts and accounting

Every query execution records:

- run ID, protocol revision, service/endpoint/version, exact query, filters,
  sort, UTC start/end, pagination and result cap;
- raw result count, pages attempted/retrieved, throttles, errors, retries,
  inaccessible results, and durable hash-addressed raw receipt;
- counts screened, included, excluded by reason, deferred, duplicated, and sent
  to second-person review;
- citation edges attempted and resolved;
- reviewer time, acquisition time, tool time, storage, and any paid-access cost;
  and
- changes from the frozen protocol, which invalidate the affected run until
  repeated under a new revision.

Retain failed queries and zero-result runs. A service cap, timeout, unavailable
API, paywall, or robots/access restriction is a gap, not a zero-result search.

## 11. Comparator and completeness reporting

Each comparator is frozen by content hash or immutable version/revision where
lawful, retrieval time, and inclusion rules. Map each comparator entry to a
source, an excluded record with reason, or an unresolved gap. Report separate
denominators for papers, standards, code/artifacts, advisories, and incidents.

Permitted coverage statements are bounded, for example: “93 of 100 records in
comparator X revision Y were mapped as of cut-off Z.” Prohibited statements at
this draft stage include “all attacks,” “all incidents,” “complete literature,”
or a safety conclusion inferred from no found record.

The master denominator remains `unfrozen`. This protocol, the schema, and the
seed do not change it.

## 12. Review, corrections, and exit gates

The protocol and corpus use separate reviewers for search design, claim
verification, and semantic validation where staffing permits. A reviewer may
not satisfy second-person review of their own extraction. Corrections append a
new manifest revision with a machine-readable reason and predecessor hash;
prior manifests remain available.

`WS6-G1` can change from `unknown` only when a protocol revision is frozen,
independently reviewed, executable, and accompanied by successful and failed
search-run receipts across every declared source class. `WS6-G2` additionally
requires the acquired corpus, lawful hashes or immutable revisions, claim-level
locators, supersession/correction graph, gap report, frozen comparator matrix,
semantic verifier, and second-person review of every critical claim.

Until those conditions are met, this directory is a draft control artifact and
must not promote an attack, safety, incident-attribution, performance, solver,
or historical-completeness claim.
