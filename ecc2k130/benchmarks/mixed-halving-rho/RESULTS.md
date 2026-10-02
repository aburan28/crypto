# Point-dependent mixed-halving rho: finite no-go

Protocol commit: `793ea7418055744802ed8f87ae4c78750f8c1ae0`.

Decision: **do not build the CUDA mixed-halving candidate.** Every candidate
whose optimistic free-dispatch primitive rate exceeds 26 B/s fails the frozen
finite-map admission rule. Point-dependent selection removes the scheduled
map's explicit phase correlation, but a high halving fraction makes the map
approach a permutation. Its indegree collisions fall and its first-repeat
constant rises faster than the primitive rate improves.

This is an exact degree-23 functional-graph result plus a primitive-rate model.
It is not a GPU measurement, a challenge-scale collision constant or a key
recovery.

## Controls and exact checks

The generated subgroup has order 2,095,853 and 45,562 classes under
Frobenius and negation. Both selectors passed 8,383,408 covariance checks with
zero failures. Sixteen deterministic random functional graphs of the same
quotient size had:

| statistic | random control |
|---|---:|
| mean exact quotient 2-cycles | 0.5625 |
| mean basin fraction ending in a cycle of length at most 8 | 0.046182 |
| median across maps of mean first repeat / sqrt(N) | 0.955896 |
| mean indegree collision pairs | 22,779.1 |

The frozen gates derived from those controls were at most 2.125 mean
two-cycles, at most 0.112363 short-cycle basin fraction, and a normalized
first-repeat median in `[0.716922, 1.433843]`.

## Rows above the raw 26 B/s line

`raw B/s` is the harmonic primitive model using the measured 28.953 B/s
halving and 20.1343255 B/s historical table-add rows. It assumes free dispatch
and compatible representations. `repeat-adjusted diagnostic` scales that raw
model by the random-control first-repeat median divided by the candidate's.
It is shown only to make the finite failure's direction visible; it is not a
challenge projection.

| selector | half threshold | actual half fraction | raw B/s | mean quotient 2-cycles | repeat / sqrt(N), median | min indegree collision pairs | repeat-adjusted diagnostic B/s | failed gates |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| canonical hash | 12/16 | 0.748431 | 26.079 | 4.1875 | 1.593427 | 9,945 | 15.645 | two-cycle, repeat |
| canonical hash | 14/16 | 0.874303 | 27.442 | 1.2500 | 2.329882 | 5,320 | 11.259 | repeat |
| canonical hash | 15/16 | 0.936943 | 28.175 | 0.2500 | 3.116542 | 2,762 | 8.642 | repeat |
| necklace7 | 12/16 | 0.742154 | 26.015 | 3.6875 | 1.538867 | 10,105 | 16.160 | two-cycle, repeat |
| necklace7 | 14/16 | 0.874303 | 27.442 | 0.8750 | 2.552194 | 5,315 | 10.278 | repeat, add balance |
| necklace7 | 15/16 | 0.941047 | 28.224 | 0.3125 | 3.641839 | 2,581 | 7.408 | repeat, add balance |

The pure-halving 16/16 control has zero indegree collision pairs and is a
permutation, as required. The pure-add 0/16 controls reproduce the additive
inverse-edge problem, with about 65 quotient two-cycles and large short-cycle
basins. Intermediate rows interpolate between those two failures; none passes.

## Interpretation

The important quantity is not how many cheap permutation steps can be counted.
The table-add branch supplies the map's merging. Making it rare lowers the
number of preimage pairs and increases the work needed to reach a repeat. The
12/16 rows are the closest to the requested raw line and already miss both the
random-like two-cycle and repeat gates. Raising the halving fraction fixes the
two-cycle count but makes the repeat constant substantially worse.

No CUDA run is admitted by the preregistered rule. The exact rows and all
per-seed graph statistics are retained in `f23-results.jsonl`;
`summarize.cpp` recomputes the accepted `summary-native.json` and the decision.
An earlier local Python-derived summary is disclosed in `NATIVE-REPLAY.md` but
is not an accepted or committed research execution path. The native replay
agrees with that superseded local derivation in every gate and aggregate field
to 12 decimal places.

## Reproduction

```sh
clang++ -O3 -std=c++17 -Wall -Wextra -Werror -Wno-unused-function \
  benchmarks/mixed-halving-rho/mixed_halving_host.cpp \
  -o build/mixed-halving-host
./build/mixed-halving-host \
  benchmarks/mixed-halving-rho/f23-results.jsonl 16
clang++ -O2 -std=c++17 -Wall -Wextra -Werror \
  benchmarks/mixed-halving-rho/summarize.cpp \
  -o build/mixed-halving-summarize
./build/mixed-halving-summarize \
  benchmarks/mixed-halving-rho/f23-results.jsonl \
  benchmarks/mixed-halving-rho/summary-native.json
```
