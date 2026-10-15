# Prime Degree Frontier

Find the largest odd prime degree for which you can produce a complete,
independently verified isogeny from each registered prime-field curve.
Each exact ICV1 source model has its own leaderboard. All registered standard prime
models are listed, including models whose construction or replay adapter is
currently unavailable. No curve is substituted or silently omitted.

The roster is the pinned catalogue's prime-field models with nonempty
`standard_names`. Names that describe the same exact model share its ICV1 row.
Registered research descendants and benchmark models do not enter the standard
curve leaderboard. The complete catalogue remains available through `isogeny curves`.

This is an open-ended computational challenge. There is no finite largest
isogeny degree in mathematics. A result here is a record for the pinned
challenge catalogue and verifier, not a worldwide record or an upper bound.

## Run the dedicated executable

The challenge is built into `isogeny`, distributed by the repository's
ReleaseMe archive workflow. It needs no Cargo, checkout, companion program,
external catalogue or rules file at runtime. Source builds are maintainer work.
Actual publication of this new command must be verified after it lands.

```sh
isogeny challenge policy > policy.json
isogeny challenge plan --curve all --from 10454 --to 65537 > ladder.json

# A short end-to-end control; scores only its own source curve.
isogeny challenge run --curve p192 --ell 13 --out fresh-control

# Attempt the recorded P192 degree before attempting a higher rung.
isogeny challenge run --curve p192 --ell 10453 --out fresh-p192-10453

# Judge a map produced by any algorithm; directory must be fresh.
isogeny challenge judge --input submission --out fresh-judgement

# JSON paths are resolved relative to this manifest's directory.
printf '%s\n' '["fresh-p192-10453/construction", "other-submission"]' > submissions.json
isogeny challenge leaderboard --input submissions.json --out fresh-board
```

`policy` pins the rules and exact embedded catalogue by SHA256. Their combined
hash identifies a challenge season. `plan` emits a descending-degree shortlist
for every requested curve. Its low-order Frobenius filter is not exhaustive
existence testing, and no candidate earns a score. Missing adapters and orders
remain visible. The default ladder starts above the previously recorded P192
degree 10453; other curves retain separate incumbents and scores.

## Admission and ranking

Submit the existing `isogeny search` artifact format: `search.json`, the
canonical curve directory, `ell-DEGREE.json`, and its stderr file. Construction
receipt hashes must match the actual bytes. Include the full monic kernel,
normalized rational-map numerator/denominator and codomain. A stored PASS,
score, ideal description, point evaluation or long chain is insufficient.

The judge snapshots bounded input files into a fresh directory and reruns the
independent verifier in a separate process. It checks the registered source,
prime degree, kernel torsion and subgroup closure, codomain, exact rational-map
identity, subgroup order and twenty public scalar transports per map. One
verified Frobenius eigenline is enough; PARTIAL/COMPLETE enumeration coverage
stays explicit. No custom curve, twist or derived walk source can replace the
registered source model. Every leaderboard build performs these checks again.

Largest verified prime degree wins **within its exact source curve**. Equal
degrees share rank. There is no speed tiebreak or cross-curve aggregate score.
Rows without accepted submissions remain UNSCORED. Timeouts, memory limits
and unavailable resource monitoring are censored and have no score. Invalid
maps are rejected. All attempt receipts remain available.

## Resources, representation and current limits

Independent replay has a maximum 3600-second budget and an 8 GiB sampled
resident-memory cap; missing monitoring fails closed. `judge --timeout N`
may request a shorter replay budget. Map JSON is limited to 256 MiB, summary
JSON to 1 MiB, stderr to 8 MiB. Resource figures are sampled process receipts,
not peak-memory proofs or isolated speed measurements. Snapshot copying and
bounded metadata handling precede the timed replay process.
The child waits for a positive RSS sample before parsing. A completion marker
ends RSS polling after all parsing and algebra, including rejected-map paths;
the process deadline continues through output and operating-system exit cleanup.

The bundled construction adapter uses at most 1800 seconds and 8 GiB per
constructor. Its current degree bound is 1000000; the judge accepts larger
u64 prime degrees when the complete explicit artifact fits its byte and replay
limits. These are implementation bounds, not nonexistence results. External
construction costs are not independently attested by the judge. They are
reported separately and do not enter the degree ranking. Full construction,
setup, replay and total costs must accompany any later efficiency claim.

Compact representations, custom curves and multiplication/composite maps are
excluded from this leaderboard. A future compact-representation track needs
its own fixed curves, representation contract and independent evaluator; it
has no accepted entries under these rules.

The immutable executable rules are in
[`tools/isogeny-cli/challenge.json`](../../../tools/isogeny-cli/challenge.json).
The prior P192 degree-10453 result is a historical anchor, not an automatically
admitted season score. It must be freshly replayed through the challenge judge.
