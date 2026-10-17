# Prime Degree Frontier

Find the largest odd prime degree for which you can produce a complete,
independently verified isogeny from each registered standard prime-field curve.
Each exact ICV1 source model has its own leaderboard. All registered standard prime
models are listed, including models whose construction or replay adapter is
currently unavailable. No curve is substituted or silently omitted.

The roster is the pinned catalogue's prime-field models with nonempty
`standard_names`. Names that describe the same exact model share its ICV1 row.
Registered research descendants and benchmark models do not enter the standard
curve leaderboard. The complete catalogue remains available through `isogeny curves`.

The [starting leaderboard](leaderboard.json) has all 54 standard prime models.
Fresh replay admits the historical P192/10453, P224/1471 and P521/7 entries;
the other rows remain unscored. The [native evidence report](../../../research/isogeny_degree_challenge_20261015/README.md)
retains the exact source/target identities, rules, raw maps and certificates.
This starting board is a historical season snapshot. The installed release
emits its own season through `challenge policy`; replay the raw submissions
again when building a leaderboard for that season. Scores do not transfer
between seasons automatically.

This is an open-ended computational challenge. There is no finite largest
isogeny degree in mathematics. A result here is a record for the pinned
challenge catalogue and verifier, not a worldwide record or an upper bound.

## Run the dedicated executable

The challenge is built into `isogeny`, distributed by the repository's
ReleaseMe archive workflow. It needs no Cargo, checkout, companion program,
external catalogue or rules file at runtime. Source builds are maintainer work.
[Release v0.1.2390](https://github.com/aburan28/crypto/releases/tag/v0.1.2390)
is published from commit `bd7eb117f094a336d3e3e8f31facb0c6b86f95bc`.
Fresh downloads of all four dedicated archives matched the release checksums;
the downloaded Apple-silicon binary passed native public-map admission controls
from outside the checkout. The [publication receipts](../../../research/isogeny_world_frontier_20261016/replay-optimization-v1/official-release-v0.1.2390/README.md)
retain the exact assets, source and platform coverage.

For Apple silicon, run this in a fresh working directory:

```sh
set -e
mkdir isogeny-2390
cd isogeny-2390
curl -fLO https://github.com/aburan28/crypto/releases/download/v0.1.2390/isogeny-v0.1.2390-aarch64-apple-darwin.tar.gz
printf '%s  %s\n' \
  cc88bac7575f51ea56dafeed98e40520cab5d5daac818d30e1bb4c293b6ec779 \
  isogeny-v0.1.2390-aarch64-apple-darwin.tar.gz | shasum -a 256 -c -
tar -xzf isogeny-v0.1.2390-aarch64-apple-darwin.tar.gz
cd isogeny-v0.1.2390-aarch64-apple-darwin
./isogeny --version
./isogeny curves
./isogeny --help
```

The version JSON reports the CLI package version `0.1.0`; its `git_commit`
identifies the exact source of rolling repository release `v0.1.2390`.
The release also provides these dedicated archives and a `SHA256SUMS` file:

| Platform | Archive |
|---|---|
| Linux x86-64, glibc | `isogeny-v0.1.2390-x86_64-unknown-linux-gnu.tar.gz` |
| Linux x86-64, musl | `isogeny-v0.1.2390-x86_64-unknown-linux-musl.tar.gz` |
| Apple silicon | `isogeny-v0.1.2390-aarch64-apple-darwin.tar.gz` |
| Windows x86-64, MSVC | `isogeny-v0.1.2390-x86_64-pc-windows-msvc.zip` |

After extracting the matching archive, run its `isogeny` executable
(`isogeny.exe` on Windows). From that directory:

```sh
./isogeny challenge policy > policy.json
./isogeny challenge plan --curve all --from 10454 --to 65537 > ladder.json

# A short end-to-end control; scores only its own source curve.
./isogeny challenge run --curve p192 --ell 13 --out fresh-control

# Attempt the recorded P192 degree before attempting a higher rung.
./isogeny challenge run --curve p192 --ell 10453 --out fresh-p192-10453

# Judge a map produced by any algorithm; directory must be fresh.
./isogeny challenge judge --input submission --out fresh-judgement

# JSON paths are resolved relative to this manifest's directory.
printf '%s\n' '["fresh-p192-10453/construction", "other-submission"]' > submissions.json
./isogeny challenge leaderboard --input submissions.json --out fresh-board
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
COMPLETE two-line coverage requires two distinct normalized kernels and two
distinct Frobenius roots. Repeating one valid map cannot establish completeness.
For a partial map with repeated roots, the judge leaves the total eigenline
count unestablished.

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
