# Protocol: stride targets for the ecbench rootfind engine (2026-10-14)

Frozen before any claim data.  Conductor T-1.  Class: **engineering** (a
constant moves; no exponent can, see §5).

## 1. Change under test

`targets=stride` (`Targets::Stride` → `TargetSource::Stride`,
`src/cryptanalysis/ic_boundary.rs`): the relation loop presents
`R₀ + [k]S`, `k = 0, 1, …`, with `R₀ = [a]G + [b]Q` and one
`S = [s]G + [t]Q` drawn per run, and the guard unchanged.  Set-up is four
scalar multiplications; every target is one addition.  The incumbent
`targets=walk` draws 16 jumps and up to 16 restart offsets at two scalar
multiplications each.  `walk` is unchanged, so every committed session
replays bit for bit.

## 2. What the exploratory runs showed (not evidence)

Scratch sessions on other target seeds (777001–777003, 2026-10-14) fixed
the arms below and are not reused:

- n = 13, 8–32 targets: stride lowered the d = 5 arm from S ≈ 25 to ≈ 13;
  among d = 3, 4, 5, 6 with the stride, d = 4 was lowest (≈ 11.8; rho-strong
  ≈ 14.5).  **d = 4 was chosen on these data**; H2 tests it on fresh seeds.
- n = 19, 8 targets: rho-strong S ≈ 3.1; every IC arm 23–51; stride and
  walk at d = 5 indistinguishable.

## 3. Sessions (specs committed with this file)

| session | spec | curve | targets × rounds |
|---|---|---|---|
| A (claim) | `docs/ecbench/specs/claim-koblitz-rootfind-stride-n13.json` | Koblitz a = 0, n = 13 (r = 2003) | 32 × 3 (+1 warm-up), alternate order |
| B (scaling check) | `docs/ecbench/specs/check-koblitz-rootfind-stride-n19.json` | Koblitz a = 0, n = 19 | 16 × 2 (+1 warm-up) |

Runner: `ecbench run --wait --cpus none --allow-busy` (L0; this container
has no CPU reservation), then `ecbench verify --replay-all`, `ecbench
table`, and `ecbench compare --a <baseline> --b <candidate> --save`
(reported ratio is candidate ÷ baseline, paired by target, with its
bootstrap 95 % interval).

## 4. Hypotheses and decision rules

- **H1 (A).** `ic-negfold-m2-d5-stride` ÷ `ic-negfold-m2-d5-walk`: the
  stride lowers S of the incumbent iff the interval's upper end is < 1.
- **H2 (A).** `ic-negfold-m2-d4-stride` ÷ `rho-strong`: complete index
  calculus is below rho-strong on this panel iff the upper end is < 1.
  rho-strong's S is a lower bound (unpriced canonicalisations), so a pass
  holds a fortiori; a fail is reported as a fail.
- **H3 (B, descriptive).** Report `…-d6-stride` ÷ `…-d6-walk` and each
  IC arm ÷ `rho-strong`.  Prediction: stride ÷ walk within 0.9–1.1 and IC ÷
  rho > 1, i.e. the stride removes a fixed set-up cost and nothing that
  grows with r.
- Every record must verify and every replay must match; the A/A control
  (`rho-strong-aa` ÷ `rho-strong`) must include 1, or the session is
  reported as unreliable.

## 5. What cannot be claimed

m = 2 decomposition costs about `|F|² + N/|F|` group operations, at best
`N^{2/3}` against rho's `N^{1/2}`; the stride changes only the constant.
No statement about larger fields, ECC2K-130, or security follows.
