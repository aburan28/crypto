# Stage 24: independent retained-mathematics replay

## Outcome and boundary

Stage 24 independently replays the mathematical witnesses retained by the
Stage-23 successor-03 terminal-evidence bundle. The replay is a separate Python
implementation of the binary field, Koblitz group law, algebraic factor base,
signed-Frobenius projection, relation reconstruction, modular elimination, and
terminal point checks. It invokes a 13-line Rust helper only for the official
`blake3` 1.8.7 hash primitive. It does not import or execute the archived
producer, retained binaries, bundled verifier, or producer verification
booleans.

The canonical replay completed 1,251 named check units with no break. Its
narrow conclusion is:

```text
retained_mathematical_witness_replay_completed = true
independent_mathematical_payload_replay_completed = false
scientific_measurement_admitted = false
external_portable_verification_satisfied = false
independent_external_reproduction_satisfied = false
full_cost_gate_passed = false
koblitz_index_calculus_sota = false
```

The broad mathematical-payload flag stays false because the archive does not
retain the hidden rho states needed to replay its exact trajectory. SAT proof
and conflict provenance are also absent. This local replay is not an external
host reproduction, novelty review, asymptotic result, cryptographic-size
measurement, or SOTA finding.

## Reproduction

Build the separately locked hash helper in a disposable directory, then pass
the helper executable and any relocated successor-03 bundle to the replay CLI:

```text
cp -R \
  research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-24-independent-math-replay-20260910/blake3-helper \
  "$TEMPORARY_BUILD_ROOT/blake3-helper"
cargo build --offline --locked --release \
  --manifest-path "$TEMPORARY_BUILD_ROOT/blake3-helper/Cargo.toml"
python3 -B scripts/koblitz_stage24_math_replay.py \
  --bundle "$STAGE23_SUCCESSOR_03_BUNDLE" \
  --blake3-helper "$TEMPORARY_BUILD_ROOT/blake3-helper/target/release/stage24-independent-blake3"
```

Run the canonical and focused tamper checks with:

```text
python3 -B scripts/test_koblitz_stage24_math_replay.py \
  --bundle "$STAGE23_SUCCESSOR_03_BUNDLE" \
  --blake3-helper "$TEMPORARY_BUILD_ROOT/blake3-helper/target/release/stage24-independent-blake3" \
  --expected-result \
    research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-24-independent-math-replay-20260910/replay-result.json
```

The CLI records identities rather than input paths, so neither the bundle nor
the helper needs to occupy its production location. The compiled helper hash is
a receipt for this host build; the helper source and locked dependency hashes
are the portable construction identity.

## Bound inputs and replay implementation

- Stage-23 successor-03 bundle seal SHA-256:
  `d2b096d39b9158e8da844d89c58732b83dfb9a38ebe08df4f1229a88a581ca5c`.
- Bundle manifest SHA-256:
  `0dc95216c59a5cd8f9aea6346c4336711f0fd4ef7b1c75481da9944118e5af48`.
- Original run seal SHA-256:
  `41a9585d975ae631d24b17963ba3de131e8ce9f309648abc65777d71673a7ada`.
- Original run inventory SHA-256:
  `de7e699d2e97f65932492a8c29f054d578b620211cda9a44a0c6749d97e01562`.
- Execution source commit:
  `54c582009762603a3311c218e4f67f0c4cf32825`.
- Execution source tree:
  `f309b301c84f3eb4026a3e71b069c9b38136e7e2`.
- Replay CLI SHA-256:
  `932c587a940454ddd963ce70335d779b7b3522659683c80ae5452987668616af`.
- Focused test SHA-256:
  `8540500cd03f15f763dd613268a8011fcc707e99add1cf6f7250521c9f11014c`.
- Canonical replay-result SHA-256:
  `f8b2c6c51b6c4420d2b73617e3dc48d8017b6a57116263bfbe93677ac2570412`.
- Helper `Cargo.toml` SHA-256:
  `549cae7bed8187e5be896384cc82393d20ab138b32159355a152e73e8781de0a`.
- Helper `Cargo.lock` SHA-256:
  `ed0b99cb0e167df661f419204b3f527760ea746d2cae6b23e80d94cc45731263`.
- Helper `src/main.rs` SHA-256:
  `6ee2ac22f677b4a758458df02cc136db33a6e873a353a7a7da04999b11a90830`.
- Canonical local helper executable SHA-256:
  `ec21c30b4d859cbcb4bbcc5d9606136ff04e72daf3bb6ff91703b30dee4b3243`.

## Check-count derivation

The 1,251 check units are fixed by the retained payload rather than wall time:

| Check class | Count |
|---|---:|
| Field, curve, factor-base, projection, target-stream, and identity checks | 18 |
| Recomputed relation-attempt targets `[a]G + [b]Q` | 437 |
| Three checks for each admitted witness: decomposition, projected row, coefficient binding | 756 |
| Target binding and relation count for five rows | 10 |
| Rank, scalar, IC point, rho point, scalar match, and rho ledger for five rows | 30 |
| **Total** | **1,251** |

Additional aggregate, hash, discovery-control, timing, resource, and claim-boundary
assertions are required by the implementation but do not inflate this stable
count.

## Field, curves, and factor base

The independent field is

```text
F_(2^23) = F_2[z] / (z^23 + z^5 + 1)
irreducible bitmask = 0x800021 = 8388641
```

For `K_0: y^2 + xy = x^3 + 1`, the replay derives:

- Group order `8,383,412 = 4 * 2,095,853`.
- Prime subgroup order `2,095,853` and cofactor `4`.
- Generator `(7,502,454, 6,195,881)`.
- Frobenius eigenvalue `lambda = 93,194`.
- `pi(G) = [93,194]G` and
  `93,194^2 + 93,194 + 2 = 0 mod 2,095,853`.

The target-independent `K_1` discovery control independently gives group order
`8,393,806`, prime subgroup `4,196,903`, cofactor `2`, generator
`(2,796,197, 1,282,171)`, and Frobenius eigenvalue `368,914`.

The factorization of `X^23 - 1` is represented by bitmasks
`[0x3, 0xae3, 0xc75]`, of degrees `[1, 11, 11]`. Selected indices `[0,2]`
give divisor `0x149f = 5279`, dimension 12, and linearized exponents
`[0,1,2,3,4,7,10,12]`.

The replay reconstructs 4,096 abscissae, 4,281 rational points, 95 signed
Frobenius orbits, 4,279 distinct cofactor-projected point keys, and 93 projected
signed-Frobenius columns. For all 4,281 points it checks

```text
[4]P = (-1)^s pi^k(R_o).
```

The factor-base identity BLAKE3 is
`47e048405e8dd0ecf8f5bae60805e601f9d18c960b081d072ea0f1f25a76f320`.
Both alternative `K_0` and both `K_1` candidates reproduce their archived
point/orbit counts and the frozen selection rule.

## Targets and row results

The target stream reproduces accepted draw counters `[1,2,3,4,7]`, candidate
attempt counts `[2,1,1,1,3]`, every packed point, full target ID, IC seed, and
rho seed. Every target is on `K_0`, nonidentity, and satisfies
`[2,095,853]Q = O`.

| Row | Target ID | `(x,y)` | IC seed | Rho seed | `d` | Attempts / Unknown | Relations / rank / nullity | Rho iterations / restarts |
|---:|---|---|---:|---:|---:|---:|---:|---:|
| 1 | `09340a19ba52d5b00e6aefce98b5a4e93cbfcc95f2774621490239e2b66fce84` | `(1885398,5393464)` | `4683119321296053311` | `15161383483049228692` | `1411444` | `79 / 36` | `43 / 43 / 51` | `1304 / 6` |
| 2 | `522066fcd091ee2b77a311267bafececcbe0ea2326d0e03ad179e1299a512eb7` | `(3229622,8259389)` | `1620815326269526617` | `16245391147071081286` | `785405` | `103 / 43` | `60 / 60 / 34` | `396 / 1` |
| 3 | `49a9376307924b4ef222dbba7750ee6617a673aafb3e69e7afcc7d3149af9a31` | `(7300722,5060417)` | `2905152658302948179` | `5342667305116464321` | `1336591` | `62 / 22` | `40 / 40 / 54` | `385 / 2` |
| 4 | `2744c699e6af05cf4abd82def616a6d54cd93ca0f231d3e5dfb41fb3a3d7b58b` | `(1884480,4110082)` | `10297604926406953080` | `12015459947305417849` | `1331121` | `120 / 52` | `68 / 68 / 26` | `494 / 3` |
| 5 | `e4a428e098bfd3b16f659d7bd081227aeacd01eddcac50cf5771d6ddc5da2767` | `(995444,2646432)` | `5031350938967679650` | `16467523183750131884` | `179888` | `73 / 32` | `41 / 41 / 53` | `168 / 1` |

Every one of the 437 attempt targets matches `[a]G + [b]Q`. All 252 admitted
two-point decompositions match their independently materialized factor-base
points. Every projected summand, sign, Frobenius exponent, coefficient, and
93-column relation row matches. Fresh modular elimination independently pins
the five target coordinates despite matrix nullities, and both the IC and rho
reported scalar satisfy `[d]G = Q` in every row.

## Aggregate accounting

The replay totals 437 solver calls, 252 admitted relations/models, 185 capped
`Unknown` outcomes, zero refuted/UNSAT outcomes, zero invalid models, and
28,422,672 reported conflicts.

Rho accounting reconstructs 2,747 iterations, 13 restarts, 442 coefficient
draws and setup scalar multiplications, 221 setup additions, 8,241 walk
additions and partition hashes, 8,254 canonicalizations, 189,842 Frobenius maps
and negations, 13 collisions, eight failed collisions, and five final scalar
multiplications. All per-row progress and ledger identities balance.

The five rho processes report 10,916,874 ns of target validation, 139,199,794 ns
of setup, 266,998,622 ns of walk time, 1,773,249 ns of candidate verification,
and 421,731,124 ns end to end. Those component values sum to 418,888,539 ns,
below the end-to-end total. This is ledger consistency, not an independent clock
trace.

Raw receipts reconstruct 14 processes, 648.3490439999999 user-seconds,
21.278418000000002 system-seconds, 669.627462 total core-seconds,
1,144.8792450850597 summed process-wall-seconds, and 1,194,147,840 bytes maximum
process RSS. The online and discovery-charged IC/rho ratios are respectively
`1808.0080514500391` and `1810.5960584685265`.

## Unavailable provenance

The 185 non-successful SAT attempts are all capped `Unknown`, not UNSAT. The
archive contains no clauses, assignments, proof logs, decision trails, or
per-conflict trace. Their reported conflict total can be summed, but the solver
history cannot be independently reproduced.

The archive also omits rho jump-table coefficients, initial coefficients,
per-step tortoise/hare states, partition selections, collision coefficient
pairs, and failed-candidate derivations. Terminal point equations, milestone
ordering, and aggregate ledgers pass; the hidden rho trajectory does not become
proved by those checks.

The five focused tamper controls alter an attempt target, projected relation
coefficient, recovered scalar, coherent modular right-hand side, and rho
walk-addition ledger. All are distinguished or rejected. These controls do not
fill the missing proof or trajectory data.
