# Frozen public points before IC execution

The protocol and build receipt were published in draft PR #1353 before these
four points were generated. Each point came from the first public
hash-to-curve/cofactor output of the frozen native `koblitz_rho_fixture`:

```sh
timeout 30 target/release/examples/koblitz_rho_fixture \
  <n> 0 signed_frobenius 1 strong <rho-seed> hash:<public-hash-seed>
```

The four fixture commands exited 0. Their raw one-line JSON and stderr are
under `targets/`; all four fixture rows state `verified: true`,
`scope: public_hash_unknown_scalar`, `target_kind:
public_hash_to_curve_cofactor`, and `published_fixture_scalar: null`.
`targets/*.jsonl` contains only the point `[x,y]` passed to IC. No IC run
occurred before this input freeze. The fixture's recovered scalar is retained
in its raw JSON as a later cross-check, but is not passed to IC.

| Phase | `n` | Hash seed | Fixed rho/rank seed | Frozen Q |
| --- | ---: | ---: | ---: | --- |
| Pilot | 41 | 41261107 | 410041 | `[1119268120096,860021860491]` |
| Held out | 41 | 41261207 | 410041 | `[1071506060992,898053054019]` |
| Pilot | 53 | 53261107 | 530053 | `[4279047133354735,1722162599174744]` |
| Held out | 53 | 53261207 | 530053 | `[7960849849661793,7443722527872608]` |

Each fixture was checked for its exact `n`, `a=0`, hash and batch seeds,
nonidentity two-coordinate point, target kind, scope and verified scalar.
The point-coordinate strings had no matches in this sparse checkout's
`experiments`, `research`, `docs`, or `examples` outside this panel. GitHub
code search for each of the four x coordinates returned no match in the
repository's indexed code at freeze time. This is a duplicate screen of the
available current corpus, not a search of every historical Git commit.
`TARGET_SHA256SUMS` pins the raw fixtures, empty stderr files and exact Q
files. The held-out Q files remain unopened by the IC producer until the
pilot selection is frozen and published.
