# ISO-1 population measurement in 192–252-bit fields

Registered before the population run, 9 October 2026. The user requested
measurement of the discovered weakness at 192–252-bit field sizes. This round
retains two distinct questions: how often an independently sampled curve passes
the proved necessary condition, and how often an explicit weak model is found
in its isogeny class. A passed conductor test is an admission, not a positive
class label. A capped or exhausted degree-2/3 search leaves the class unresolved.

## Fixed experiment

- Fields are `K=F_(p^(2n))`, with base `F_(p²)` and odd degree `n`.
- Six degree-6 fields have exact integer bit lengths 192, 204, 216, 228, 240,
  and 252. Their primes are respectively 4294967291, 17179869143, 68719476731,
  274877906899, 1099511627689, and 4398046511093: the largest primes below
  `2^32,2^34,...,2^42`.
- Degree-10 endpoint fields use `(p,bits)=(602233,192),(38543917,252)`;
  degree-14 endpoint fields use `(13421,192),(262139,252)`. Each prime is the
  largest prime below the integer `2n`th root of `2^bits-1`.
- Draw 64 curves per degree-6 field and 32 per degree-10/14 field: **512
  independent curve fixtures**. Sampling is uniform in ordered pairs
  `(u,v)` of distinct nonzero elements of `K`, giving
  `E:y²=x(x-u)(x-v)`. This samples full rational 2-torsion models;
  it is not uniform in traces or isogeny classes. Repeated traces stay rows.
  Sampling is independent of the norm-one construction and conductor outcome.
- Every field initializes PARI's field generator from seed `2026100900`.
  Curve seeds are `202610090000 + 10000*field_number + sample_number`.
  All exact moduli and coefficient vectors are recorded.
- Ordinary status is `p∤t`; the necessary admission is
  `t ≡ ±(Q+1) mod 16`, equivalent here to `4|f_pi`. Record the exact
  2-adic depth of the Frobenius-order conductor without factoring its odd part.
- Independently count each source using PARI/GP `ellcard`. Verify Hasse,
  full rational 2-torsion, and `[N_E]P=O` on three generated points.
  Record direct norm weakness, full rational 4-torsion, and a rational
  2-isogeny neighbor with full 4-torsion after choosing the admitted trace sign.
- For each ordinary admitted class, breadth-first search rational degree-2
  and degree-3 isogenies, retaining full-2-torsion vertices. At most **256
  distinct j-invariants** are tested; the search has a **90-second** cap.
  A witness must pass the direct norm criterion and have a verified route.
  Search exhaustion of this restricted component is not a certified zero.
- Every retained isogeny is evaluated on three source points; all target
  2-torsion roots are checked against its equation. A found witness is freshly
  point-counted and must have the selected source trace.
- One separately labeled forced-positive control per field starts at the
  degree-2 neighbor of a norm-one weak model. These check that the witness
  search recognizes a known weak class; they are excluded from population rates.
- A small-field `p=7,n=3` validation compares witness labels with the complete
  exact census, and independently counts every retained route vertex.
- Native Rust orchestration, four GP workers, 256 MiB initial GP stack,
  240-second parent cap per sample. All statuses, caps, stderr, raw output,
  source and binary hashes are retained. The process cap is separate from
  GP's search cap. Wall times are resource receipts on a shared host, not
  controlled performance comparisons.
- Report Wilson 95% intervals for admission proportions. For admitted classes,
  if `w` have verified weak witnesses, `z` certified zeros, and `u` unresolved,
  the identified precision interval is `[w/(w+z+u),(w+u)/(w+z+u)]`.
  Do not replace unresolved labels with false positives or transfer the
  small-prime census precision to these fields.

## Requirement status gates

| Requirement | Completion evidence |
| --- | --- |
| 192–252-bit field sizes | Exact cardinalities, bit lengths, irreducible moduli |
| Independently sampled population | 512 frozen, non-selected source fixtures and count receipts |
| Weakness measurement | Admission, torsion, direct norm test, and bounded witness yield recorded separately |
| Class prediction precision | Exact labels or explicit unresolved bounds; remains partial if admitted labels are unresolved |
| Larger extension degrees | Degree-10 and degree-14 endpoint measurements |
| Reproducible publication | Native runner/checker, source freeze, report, diagram, plot, PDF, same PR |

The previous complete small-prime census and its ongoing worker remain separate
evidence. Their ordinary-row population and probabilistic point-count labels
must not be pooled with this independently sampled curve population.
