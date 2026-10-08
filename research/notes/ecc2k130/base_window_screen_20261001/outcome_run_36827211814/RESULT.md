# n=41 four-window screen: no window opportunity

**Decision: `NO_WINDOW_OPPORTUNITY`.** Class: diagnostic. This is not a selector, not a method speedup, and not an ECC2K-130 result. Common operation-equivalent `S`, a second host, and n=131 transfer stay unset.

The one-shot hosted run is Actions `36827211814`, `workflow_dispatch` on main `4978478120ae83a322c221e59d65a358d5510f29`. `validate`, `measure`, and `replay` all succeeded. No second dispatch was started. The input freeze and source lock are the ones already on that commit. Curve `icv1-f2m41-tm2308219-7f48b14a`, K=255, L=1,024, five blocks.

## What was charged

Forty-five children were planned and recorded. The independent replay verified 40,960 compact target logs, 5,120 rho target logs, and 40 rank replays. Timing was eligible: the machine was uncontended, and every window's A/A interval contains 1. Each compact ratio below is generator-plus-compact child CPU divided by the same block's rho child CPU. A value above 1 means that arm used more CPU than rho.

| Window | Geometric mean | 95% interval | A/A geometric mean |
| ---: | ---: | --- | ---: |
| 0 | 1.3797151711313378 | [1.338080223061073, 1.4226456087178052] | 0.9939879383861577 |
| 1 | 1.3963992292185985 | [1.3387852330904553, 1.4564926167142365] | 1.0202685843542891 |
| 2 | 1.3987290065775213 | [1.3629045689237829, 1.435495102482666] | 1.0249016303604515 |
| 3 | 1.423787713341646 | [1.3626974448703122, 1.4876166828474233] | 1.0053895323699447 |

No fixed window has an interval below 1, so none is a crossover candidate.

The post-hoc per-block minimum over the four windows has geometric mean 1.3705041917887892 and 95% interval [1.334355540712842, 1.4076321358153336]. Every block's minimum is above 1, and the interval's lower end is above 1. That minimum pays each measured policy and pays zero selection cost. It is an optimistic lower envelope, not a policy that can be deployed. Because even that envelope stays above rho, a nonnegative-cost selector restricted to these four bases cannot cross on these blocks.

## Scope

Replay status is `PASS`. Screen-run SHA-256 `05d892a144d3eee6b1d7a0eed5280f2f378ab2c55bfa8ab0444069db65b2d22e`. Isolation SHA-256 `9e32da87edf6205aba31e92d68c0ff6049bf1e142d98bceb3e69bb46f5b23ce4`. The raw and replay zip hashes, and SHA-256 of all 421 replay-download members, are in `ARTIFACTS.json` and `MANIFEST.sha256`. The trace bodies stay in the Actions artifacts.

This stops the four-window family under the frozen rule. It does not price another base, a shared scan, a degree-263 descendant, or n=131. The protocol's named successor is natural high-arity descendant-native PDP yield and a non-eager index or query policy.
