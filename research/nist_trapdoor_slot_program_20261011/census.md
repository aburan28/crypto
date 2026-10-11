# NIST prime-curve trapdoor-slot census

Executed by `examples/nist_slot_census.rs`.  Every value is recomputed from `CurveParams` and the published seeds; nothing is pinned by hand.

## Group invariants (slots S2, S5, S6)

| curve | p bits | trace sign | n prime | k ≤ bound | \|D\| bits | square part of \|D\| from small factors (bits) | \|D\| residue (bits) | residue square? | residue prime? |
|---|---:|---|---|---|---:|---:|---:|---|---|
| P-192 | 192 | + | true | > 1000 | 194 | 1 | 184 | false | true |
| P-224 | 224 | + | true | > 1000 | 226 | 4 | 182 | false | false |
| P-256 | 256 | + | true | > 1000 | 258 | 1 | 255 | false | false |
| P-384 | 384 | + | true | > 1000 | 386 | 1 | 386 | false | false |
| P-521 | 521 | + | true | > 1000 | 523 | 1 | 521 | false | false |

## X9.62 root-sign freedom and twist (slots S3, S5)

| curve | seed → b holds | p mod 4 | sign choice = twist | published root smaller | twist order prime | twist leak bits (< trial bound) | twist residue bits | residue prime |
|---|---|---:|---|---|---|---:|---:|---|
| P-192 | true | 3 | true | true | false | 4.5 | 188 | false |
| P-224 | true | 1 | false | false | false | 12.2 | 212 | false |
| P-256 | true | 3 | true | true | false | 15.1 | 241 | true |
| P-384 | true | 3 | true | false | true | 0.0 | 385 | true |
| P-521 | true | 3 | true | true | false | 5.1 | 516 | false |

## Generator provenance (slot S4)

| curve | G on curve | smallest valid x | G.x is smallest valid x | G.y is smaller root | x scanned |
|---|---|---|---|---|---:|
| P-192 | true | 0x0 | false | true | 1 |
| P-224 | true | 0x3 | false | false | 4 |
| P-256 | true | 0x0 | false | true | 1 |
| P-384 | true | 0x0 | false | true | 1 |
| P-521 | true | 0x0 | false | false | 1 |

## Twist order factorisation (trial division + bounded Pollard rho)

| curve | factors (hex, exponent, probable prime) | largest proved prime factor (bits) |
|---|---|---:|
| P-192 | 0x17^1 (prime) · 0xb21642c8590b21642c8590b25d4ea13d29723c09f184109 (unsplit residue, composite) | 5 |
| P-224 | 0x3^2 (prime) · 0xb^1 (prime) · 0x2f^1 (prime) · 0x2e0273^1 (prime) · 0x268160f^1 (prime) · 0x20905fa28a700b3b10f0d8d3df7028d8feb21bc6ff (unsplit residue, composite) | 26 |
| P-256 | 0x3^1 (prime) · 0x5^1 (prime) · 0xd^1 (prime) · 0xb3^1 (prime) · 0x1e0a75640070a738557cc30f68bd56eaea65c94f98411d17ac4e16ece1a47^1 (prime) | 241 |
| P-384 | 0x1000000000000000000000000000000000000000000000000389cb27e0bc8d21ea7e5f24bb74f58851313e697333ad68d^1 (prime) | 385 |
| P-521 | 0x5^1 (prime) · 0x7^1 (prime) · 0x4277ffb^1 (prime) · 0x25e6d2cd^1 (prime) · 0x17c8b8fa594c0fc63a5c0043ab498c1762d92f18fdfe2fea8f074695615d886d81bf930a0ac77d01bf9dd8c1a1ae121dab4e860c5dc18e265de3^1 (prime) | 461 |

## |D| small factors

- P-192: [(5, 1), (11, 1), (31, 1)]
- P-224: [(3, 3), (29, 1), (79, 1), (7523, 1), (40927, 1)]
- P-256: [(3, 1), (5, 1)]
- P-384: []
- P-521: [(5, 1)]

## n ± 1 smoothness (slot S6)

| curve | n−1 small factors | n−1 residue bits | residue prime | n+1 small factors | n+1 residue bits | residue prime |
|---|---|---:|---|---|---:|---|
| P-192 | [(2, 4), (5, 1), (2389, 1)] | 175 | false | [(2, 1), (3, 1), (41, 1)] | 185 | false |
| P-224 | [(2, 2), (3, 6), (5, 1), (17, 1), (2153, 1)] | 196 | true | [(2, 1), (7, 1), (19, 1)] | 216 | false |
| P-256 | [(2, 4), (3, 1), (71, 1), (131, 1), (373, 1), (3407, 1), (17449, 1), (38189, 1)] | 188 | false | [(2, 1), (5, 1), (1879, 1)] | 242 | false |
| P-384 | [(2, 1), (3, 2), (7, 2), (13, 1)] | 371 | false | [(2, 2), (83, 1), (181, 1), (421, 1)] | 360 | false |
| P-521 | [(2, 3), (7, 1), (11, 1), (1283, 1)] | 502 | false | [(2, 1), (3, 2), (5, 2)] | 513 | false |

## Seed relations (slot S7)

- P-192: 84 candidates checked, 0 matched
- P-224: 84 candidates checked, 0 matched
- P-256: 84 candidates checked, 0 matched
- P-384: 84 candidates checked, 0 matched
- P-521: 84 candidates checked, 0 matched
