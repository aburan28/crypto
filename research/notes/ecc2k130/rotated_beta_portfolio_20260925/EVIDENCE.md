# Exact five-base portfolio evidence

The files here are the complete compact output of the frozen two-child
analysis at `/private/tmp/rotated-beta-portfolio-run-20260925`; the committed
copy is the durable evidence. The source and prior full-census archives are
identified by `FROZEN.json`. `outcome.json` contains all exact membership
patterns, all 31 subset unions, the fixed-order prefixes and the support-only
decision. `verification.json` is the independent all-q group-index result;
`receipt.json` contains commands, UTC intervals, exit codes, high-water RSS,
source/input and output SHA-256. The four stdout/stderr files are retained even
though they are empty.

Run from the repository root:

```sh
python3 research/notes/ecc2k130/rotated_beta_portfolio_20260925/ci_replay.py --evidence research/notes/ecc2k130/rotated_beta_portfolio_20260925/evidence
```

The replay reads only committed input archives and reconstructs the full
all-q comparison; it does not trust the local temporary directory.
