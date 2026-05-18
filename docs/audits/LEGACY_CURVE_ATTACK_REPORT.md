# Legacy Curve Attack Demonstrations

Generated with:

```sh
RUSTFLAGS=-Awarnings cargo run --quiet --example legacy_curve_attack_demo
```

This report demonstrates concrete secret recovery on the newly added legacy low-security curves. Each row uses the real curve parameters and generator, plants a deterministic reduced-size scalar, computes Q = d*G, and recovers d with bounded baby-step giant-step.

The result is a successful attack on implementations that use low-entropy or otherwise bounded private scalars on these curves. It does not pretend that a full 2^56, 2^64, 2^76, 2^80, or 2^96 generic ECDLP completes in a normal test run; the full-width estimate column gives the expected Pollard-rho scale.

| curve/group | family | secret | recovered | bound | BSGS additions | time | full-width generic rho |
|---|---|---:|---:|---:|---:|---:|---|
| secp112r1 | prime-field | 4660 | 4660 | 2^16 | 275 | 165908 us | approx 2^56 group ops |
| secp112r2 | prime-field | 9029 | 9029 | 2^16 | 292 | 161982 us | approx 2^55 group ops |
| secp128r1 | prime-field | 13398 | 13398 | 2^16 | 309 | 190367 us | approx 2^64 group ops |
| secp128r2 | prime-field | 17767 | 17767 | 2^16 | 326 | 244076 us | approx 2^63 group ops |
| secp160k1 | prime-field | 22136 | 22136 | 2^16 | 343 | 412262 us | approx 2^80 group ops |
| secp160r1 | prime-field | 26505 | 26505 | 2^16 | 360 | 418192 us | approx 2^80 group ops |
| secp160r2 | prime-field | 30874 | 30874 | 2^16 | 377 | 426215 us | approx 2^80 group ops |
| P-192 | prime-field | 35243 | 35243 | 2^16 | 394 | 439628 us | approx 2^96 group ops |
| sect113r1 | binary-field | 4951 | 4951 | 2^16 | 276 | 291842 us | approx 2^56 group ops |
| sect113r2 | binary-field | 9320 | 9320 | 2^16 | 293 | 318052 us | approx 2^56 group ops |
| sect131r1 | binary-field | 13689 | 13689 | 2^16 | 310 | 517747 us | approx 2^65 group ops |
| sect131r2 | binary-field | 18058 | 18058 | 2^16 | 327 | 555647 us | approx 2^65 group ops |
| ike-oakley-group3 | binary-field | 22427 | 22427 | 2^16 | 344 | 761758 us | approx 2^78 group ops |
| sect163k1 | binary-field | 26796 | 26796 | 2^16 | 361 | 907658 us | approx 2^81 group ops |
| sect163r1 | binary-field | 31165 | 31165 | 2^16 | 378 | 916984 us | approx 2^81 group ops |
| sect163r2 | binary-field | 35534 | 35534 | 2^16 | 395 | 964461 us | approx 2^81 group ops |

## Caveats

- The BSGS runs are exact recoveries for the planted bounded scalars shown above.
- Full-width uniformly random scalars on these curves remain generic-DLP problems at the listed Pollard-rho scale; that scale is too large for CI.
- Cofactor-bearing curves still require subgroup validation in protocols. This demo focuses on bounded-secret recovery because it is deterministic and applies uniformly across the added curve implementations.
