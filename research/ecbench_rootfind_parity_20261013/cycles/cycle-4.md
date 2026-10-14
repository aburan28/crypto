# Cycle 4 — Frobenius orbit / subtract (mostly dead ends)

Session: `sessions/cycle4` (`ECBS1h6718277fab42`), 30/48 verified
(18 errors), audit OK on verified arms.

## Hypothesis

`koblitz-orbit:divisor=0;1` + `mitm-frobenius` (and counted) should cut
cold set-up below negation-folded PairTable; `subtract` has zero
`oracle_setup`.

## Outcome

| arm | result |
|---|---|
| `ic-frob-orbit-*` | **error**: `no invariant subspace for divisor [0, 1] on this curve` (E₀ / n=13) |
| `ic-subtract-m2-d5` | **dead end**: mean S **119.8** (setup 0, relations explode) |
| `ic-mitm-neg-m2-d5` | ties negfold (same PairTable path): mean S **25.0** |
| `ic-negfold-m2-d5` | mean S **25.0** |
| `ic-f4-m3-d5` | mean S **16.5** |

## Follow-on probe (not a full cycle)

`compact-orbit-scan` is Frobenius-closed on E₀. Probe
`/tmp/probe-compact`: columns=4 + `mitm-frobenius:m=2` mean S **~19.8**
vs negfold **~22.7** vs F4 **~15.3**. Promoted to cycle-4b / cycle-5.
