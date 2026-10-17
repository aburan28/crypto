`ratios_to_rho.tsv`: mean priced S of each IC arm over mean priced S of the
cheapest constant-memory generic arm on the same curve (`rho-frob` on Koblitz
curves, `rho-neg` on prime curves) at w = 0, 0.1, 1 GAE per native operation.
Produced with `jq` from `band-*.json` (ecbench band output): a quotient of two
recorded means, nothing else. BSGS is excluded from the denominator because it
holds √r memory (PROTOCOL.md A5).
