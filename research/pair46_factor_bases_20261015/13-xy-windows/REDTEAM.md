# 13-xy-windows: red-team review

Reviewer: independent red-team subagent, 2026-10-15. No other file in this directory was changed.

## Verdict: claim upheld (null)

## Checks
- **Protocol before data.** PROTOCOL.md mtime 16:41:59, sha256 matches PROTOCOL.sha256 (f7224755...f876);
  src/main.rs 16:43, results.json 16:50. Input csv sha256 88c26329...7c12 matches the protocol.
  Arms, K, T=40, seeds, bootstrap and decision rule in the code match the protocol.
- **Reproduction.** Rebuilt and re-ran `xy_windows . 40 100 100 2000 1500` on `taskset -c 3` into a scratch
  directory. results.json is identical to the committed one in every field except wall times.
- **Same instrument across arms.** All 11 arms share one target set, one union x-mask and one
  decomposition loop (R -+ P, unordered pairs via x2 >= x1); the scan arms share one predicate function.
  Arms differ only by membership.
- **Numbers.** Every ratio, interval, scan price and control value in the headline and RESULTS.md
  matches results.json / posthoc.json (checked pooled, F and N2 groups; all joint/Xbits4/R intervals
  lie in [0.9562, 1.0257] pooled and in [0.92, 1.07] per group; X2 holds on 200/200 curves).
- **Math.** X1 (Semaev relations depend on x only), X2 (exactly one of y, p-y is even / below p/2 for
  y != 0), X3 (vanishing polynomial degree K; pole order >= K from deg div = 0) and X4 (2/q under
  equidistribution, liftable fraction ~1/2) are correct. The post-hoc binomial null 1 - q with
  q = 16/B1 is correct, and the planted index 1.23 at 2^8 matches the expected (15 + 5.76)/16.8 = 1.236.

## Caveats (do not refute)
- The positive control PC_logAP shows the detector sees a structural *deficit*; there is no planted
  control for a yield *gain*. That gain is the "advance" direction. At m = 2 and matched size, PC still shows
  the instrument is not blind to structured bases.
- The E3 rule as preregistered fails at B1 = 2^8 for all six curves. RESULTS.md discloses this. The pass is
  against a post-hoc null. The disclosure is honest, and the null is exact, not tuned. The planted threshold
  1.2 was also set under the Poisson null; at 2^8 it passes with a small margin.
- The headline "Pair-46 box counts look like those of random curves" holds only for x < 2^23, three box
  shapes and lambda ~ 16. It is not a statement about the whole curve. Corner-density exclusions
  (floor k=4..9,11; random2 k=1,2,7..9; crater k=16) are uncorrected and correlated across nested k.
- The toy yield result is m = 2 and 20-30 bit only. The pair-46 part measures equidistribution and price,
  not yield.
