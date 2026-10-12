# Fresh call-only F5 ceilings leave room for graded reuse; no cache speedup is measured

The native Linux x86-64 AVX2 screen completed and replayed **48 discovery
and 48 disjoint holdout cells, with 9,408 timed F5 calls in each split**.
Every registered n=24, batch-32 group has a 95% bootstrap **upper** bound
for the median optimistic call-only ceiling above 2.0: **4/4 discovery and
4/4 holdout groups pass the feasibility gate**. The worker and verifier
source bytes are identical across splits. This permits a separate actual
graded batch-cache implementation and timing experiment. It does **not**
measure a candidate/reference speedup.

| Split | Seed and affine family | Median `U_call` | 95% upper bound | A/A 97.5% floor | >2 upper gate |
|---|---|---:|---:|---:|---|
| Discovery | 20261011, independent | 4.692 | 4.714 | 1.008 | PASS |
| Discovery | 20261011, walk | 4.656 | 4.657 | 1.004 | PASS |
| Discovery | 3141667, independent | 4.728 | 4.742 | 1.007 | PASS |
| Discovery | 3141667, walk | 4.774 | 4.777 | 1.003 | PASS |
| Holdout | 20261018, independent | 4.398 | 4.435 | 1.006 | PASS |
| Holdout | 20261018, walk | 4.462 | 4.473 | 1.007 | PASS |
| Holdout | 4242509, independent | 4.452 | 4.459 | 1.011 | PASS |
| Holdout | 4242509, walk | 4.399 | 4.416 | 1.004 | PASS |

Here `U_call=(C+D)/(C+D-build-reduce)`. `C` is the complete inherited F5
call with its returned `Vec<F2BoolPoly>` materialized; `D` is separately
timed destruction of that returned object. Every exclusive internal phase
fits inside `C`. Exact SHA-256 output validation and small F4 reference work
are outside both timed intervals but visible in whole-process receipts.
The formula optimistically makes **all** matrix building and reduction free,
so a real cache can only do worse under the same fixed output API unless it
also changes other phases. A later candidate will add setup, fallback and
memory work. The previous validation-inclusive study's ~1.61x ceiling is
not contradicted: it timed a different boundary that charged digest work
inside `T`; neither raw run is divided into the other as a speed ratio.

Across the 56 n=24 batch observations in each split, build plus reduction
account for **78.78%** of 142,175,332,640 timed discovery nanoseconds and
**77.39%** of 144,207,697,150 timed holdout nanoseconds. The independently
recorded digest validation took 172,800,199,158 and 168,209,539,503
untimed nanoseconds, respectively. It is substantial practical work and is
not hidden, but it is not an F5 API call cost in this protocol. The phase
shares and upper bounds concern these four public synthetic quadratic cores
at n=24, not a population-wide exponent or a natural Semaev-system trace.

The inherited full-column guard used direct packed construction on 8,512
discovery calls and 8,484 holdout calls; exact sorted-row fallbacks on 896
and 924 calls remained in the fixed grid. All returned F5 report counters,
row digests, direct-unpack flags and small F4/F5 row spaces passed. The
source-defined fallback is part of the named route, not a censored or
favourably selected case.

The [discovery bundle](qualified_discovery_01/manifest.json) has manifest
SHA-256 `25ad106c29ffbe06674224ca4739e3d8cba44a84ecbeb7c4821256fc284e6557`.
The [holdout bundle](qualified_holdout_01/manifest.json) has manifest
SHA-256 `0c70e1fb0a9aa1191703ed3b156e64821b753eb4ce77bd80bb6e36b559951760`.
All 27 discovery and 30 holdout member hashes and native replays passed
after copying; the holdout binding names the passing discovery. The
identical worker SHA-256 is
`9b045d0bc9404afa2ce03c1ab98ac64dfa1793bb896cdc8dc5bd3b8b801173da`;
the verifier SHA-256 is
`5658e8707c819f669a063c2f56405f055fbf0374dafb5126e92fa4758e05a6b0`.
Both receipts were uncontended on the two SMT siblings of one physical core.
Peak whole-worker RSS was 126,056 KiB in discovery and 125,576 KiB in
holdout, including fixtures, reference calls and output validation; these
are not candidate-specific memory costs.

Nine optimized Rust tests pass, including exact F4/F5 controls, call-phase
accounting, returned digest checks, deterministic bootstrap, float
round trips, SMT topology and resource refusal. The producer and verifier
are native Rust; only the repository's existing isolation controller is
Python. No graded F5 batch cache, measured complete-call speedup, natural
relation-yield improvement, independent relation rank, calibrated IC cost,
rho crossover, curve input or key-related result is established.

**Next experiment:** implement the exact high-block cache under a new
same-binary complete-F5 batch protocol. Charge criterion recomputation on
every assignment, cached high-block setup, changing low-block updates,
output unpacking and destruction, support-signature fallbacks and peak
memory. The four tested cores establish only that this is worth measuring.
