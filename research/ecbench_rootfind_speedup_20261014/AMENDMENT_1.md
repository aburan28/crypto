# Amendment 1: `stride-renew` (frozen before its data)

Written after sessions A and B of [`PROTOCOL.md`](PROTOCOL.md) were run and
audited; those sessions and their hypotheses H1–H3 stand as registered.

## Why

Session A's d = 5 stride needed 201 trials against the walk's 119 at the
same hit rate, and finished with 14.9 rows at rank 12.1; the walk finished
with rows = rank (B: 29.6 rows at rank 23.4 against rows = rank).  Under one
fixed stride, targets are `R₀ + kS`; two pairs of relations the same number
of steps apart whose factor-base differences coincide make a four-row
dependency with `Σc = Σ c·k = 0`, which cancels the target and pins nothing.
A diagnostic with fresh random targets (seed 777010) finished with rows =
rank, as the walk does, so the stride was the outlier.

## Change

`targets=stride-renew` (`TargetSource::StrideRenew`): the stride, plus
`S ← S + T` after every target that yields a relation (`T` drawn once: two
scalar multiplications; one addition a relation), and a repeat is stepped
past instead of restarting (two progressions with different strides share at
most one point).  `stride` and `walk` are unchanged.

## Exploratory data that chose the arms (not evidence)

Seeds 777011 (n = 13, 32 targets) and 777012 (n = 19, 16 targets), one round:
`stride-renew` S = 12.74 (d = 5) and 10.98 (d = 4) at n = 13 against walk
23.63; 14.61 at n = 19 (d = 6) against walk 18.25 and rho-strong 3.13; rows =
rank in every renew arm.

## Sessions

| session | spec | targets × rounds |
|---|---|---|
| C | `docs/ecbench/specs/claim-koblitz-rootfind-stride-renew-n13.json` | 32 × 3 (+1 warm-up) |
| D | `docs/ecbench/specs/claim-koblitz-rootfind-stride-renew-n19.json` | 16 × 3 (+1 warm-up) |

Same runner, audit and comparison commands as `PROTOCOL.md` §3.  Before C
and D, the binary they use re-audits sessions A and B with `verify
--replay-all` (the `stride` path must replay unchanged).

## Hypotheses

- **H4 (C).** `ic-negfold-m2-d5-renew` ÷ `ic-negfold-m2-d5-walk`: upper
  end < 1.
- **H5 (D).** `ic-negfold-m2-d6-renew` ÷ `ic-negfold-m2-d6-walk`: upper
  end < 1.  This is the test that the gain is not only a fixed set-up
  cost.
- **H6 (C).** `ic-negfold-m2-d4-renew` ÷ `rho-strong`: upper end < 1
  (rho-strong is a lower bound; see `PROTOCOL.md` H2).
- Descriptive: renew ÷ stride in C and D; every IC arm ÷ rho-strong in D;
  rho-frob beside every n = 13 row.
- The A/A control must include 1; every record verified, every replay
  matched.

`PROTOCOL.md` §5 applies unchanged: the exponent of m = 2 decomposition
does not move.
