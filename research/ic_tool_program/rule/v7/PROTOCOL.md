# The rule's comparison at baseline v7

**Declared 2026-10-16, in the pull request that accepts v7 (R11), before
any run below.** IC_TOOL_PROGRAM.md asks for the rule's comparison at
every new baseline. v7 is v6 plus R11's streamed pair-table build.

**Everything rule v6 declared holds here, amendment 1 included, except the
rows under "What changes"** ([`../v6/PROTOCOL.md`](../v6/PROTOCOL.md)).
That covers:
- §23's six sizes, 64 targets and four A/A repeats a size, 408 rows;
- the arms, retries and stop rule;
- the claims, now with the ledger's newer `vs_rho` fields;
- the pin against §23's `T01-R1`, and each row's logarithms checked
  against §23's;
- the reference check against the strong rho, with the fixture's size gate
  widened to `2 <= n <= 63`;
- the pricing statement under AGENTS.md §5a: both arms are timed in one
  process, so no native work is unpriced on either arm, and the 0, 0.1 and
  1 GAE band gives one figure;
- the accounting, class and non-claims.

## Question

How far is rho still ahead at v7? At v6 rho was 3.78–12.06× faster cold,
and 7.90×, 7.40× and 12.06× at the three largest sizes.

## What changes from rule v6

- **The binary** is v7's, the one R11 measured: SHA-256
  `85116c2a0da19e5e2e31d7ec1ff8d58a28406bbf9ffe30d330476a67a182315a`,
  built from local commit `c2af9b7ccd92469830048dc2deaf97e1d194b27a`
  (v6 plus `R11-build-streams/candidate-v6.patch`, tree `bdd975b7`).
- **The strong rho is v7's own.** `koblitz_rho_fixture` is built from local
  commit `a475cde33904207fea38d850d7279cc87e0d27b8`: v7's tree plus the
  one-line size gate of rule v6's amendment 1. Its SHA-256 is
  `feba157ca2435b12e152126068ce2cf26c469d166b83678c891ca614df254337`.
- **What each size is read against:** §23 under the key `s23_then_v7`. The
  README sets rule v6's figure beside each one.
- **The host** is the Emerald Rapids of rule v6, kernel build
  `6.18.44-fc-v139`.

## Prediction

Not drawn from a run. If the rho arm costs what it did at v6, the cold
ratio falls by about what R11 gained in cold time: 1.134×, 1.124× and
1.071× at the three largest sizes. That gives about 7.0×, 6.6× and 11.3×.
At the smaller sizes R11 measured 1.01–1.14×.

## Success and stop

The rule makes no acceptance decision. It is a measurement.

**Stopped** if:
- a row fails its check;
- a logarithm differs from §23's;
- the pin fails;
- the reference check finds the comparison's rho slower a step than the
  strong walk at any size. That size is reported as not run.

## Commands

    icprog rule manifest --comparison v7 --ic <v7 ic> --ic-commit c2af9b7ccd92469830048dc2deaf97e1d194b27a --isolate <isolated_bench>
    icprog rule all --comparison v7 --ic <v7 ic> --ic-commit c2af9b7ccd92469830048dc2deaf97e1d194b27a --isolate <isolated_bench>
    icprog rule reference --comparison v7 --fixture <koblitz_rho_fixture> --fixture-commit a475cde33904207fea38d850d7279cc87e0d27b8 --isolate <isolated_bench>
    icprog rule claims --comparison v7
    icprog rule analyse --comparison v7 > research/ic_tool_program/rule/v7/analysis.json

These run through [`run.sh`](run.sh).
