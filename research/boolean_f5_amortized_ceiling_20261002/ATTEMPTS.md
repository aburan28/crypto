# Complete-F5 ceiling attempt ledger

## Failed discovery 01: SMT reservation refused before worker launch

The first GitHub Linux x86 dispatch used source commit
`d386b1069320e6f21ccc3646cc94cdd785a661bb` and the original frozen
protocol. The correctness job passed. The qualified job built and tested the
native worker, but `tools/isolated_bench.py run --cpus 1` refused before
starting it: **“CPU 1 shares a core with [0]; reserve the siblings too.”**
The sealed `raw.jsonl` is empty (SHA-256
`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`),
`failure.json` has `performance_admitted=false`, and no registered F5 call,
phase sample, gate or holdout was executed. This is a resource-preflight
refusal, not a mathematical negative or a timeout.

The original 24-file bundle is retained unmodified in `failed_discovery_01`.
Its 23 hashed members match manifest SHA-256
`c9a39ec563eccfb6f93649ea57a8c6d43d10820488abce09f504ac8db0890ac6`.
[Workflow 37003720117](https://github.com/aburan28/crypto/actions/runs/37003720117)
uploaded artifact ID `11225258830` with archive digest
`sha256:61959d87cd424b044cb7b80c5657acd9efbc3e36284ae7d8fbf91ecc38b7b1ce`.

The additive pre-measurement correction in [AMENDMENT.md](AMENDMENT.md)
reserves the entire physical core by reading logical CPU 2's SMT-sibling
list from sysfs. The verifier checks that the resource receipt names exactly
that recorded set. The worker remains one-threaded. The F5 algorithm,
conditional ceiling, four primary groups, 2x stop rule, runtime/evidence
caps and both seed sets are unchanged. A fresh qualified discovery is needed.
