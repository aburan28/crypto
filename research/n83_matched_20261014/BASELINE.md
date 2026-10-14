# Source and accounting boundary

Work starts from commit12321a1855d0afada849d601559d1b5eaa6e0b39 on
`codex/n83-rho-speed-20261012`, whose open result PR is
[1828](https://github.com/aburan28/crypto/pull/1828). The branch will be extended,
not replaced. New source, public fixtures, human/machine protocols and native
tests must be committed before timed solver execution. The resulting freeze
commit and executable SHA256 hashes will be recorded separately.

IC and rho are rebuilt for this experiment from the same source freeze and
release compiler profile. The preceding Oct12 timings are not measurements
of this new build and will not be combined arithmetically with it. The old
139-launch study and its comprehensive70-page report remain unchanged.

Cold starts mean fresh processes given the materialized historical K600,
without cached in-process log/index/DP state. They do not establish the cost
of discovering or constructing that historical factor base from a curve-only
input, nor OS-cold file caches. Stock IC target_ms omits input/Q-validation/
output intervals included in rho's online clock; the whole-process primary
comparison charges them for both methods. External C++ research replay is
reported separately from built-in verification already inside solver cost.

The current wide IC producer exports per-target scalar/decomposition
witnesses, rank/count/hash summaries and phase timings, but not its full600
rank-equation/factor-log transcript. Independent answer/witness replay does
not silently become independent verification of that unexported transcript.
The exact legacy fixture also lacks applicable sealed ecbench common-cost
admission and independent-host evidence. Those obligations stay open.
