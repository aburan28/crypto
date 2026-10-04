# Two-word F6 geometric closure on the n=83 K_1 field

The [frozen stage protocol](PROTOCOL.md) ran once in a release build on an
arm64 macOS 26.6 host with Rust 1.93.1. The tested curve is K_1 over the
registered n=83 field, subgroup order 8,569,786,107,849,059. It is distinct
from the K_0 n=83 confidence-gate curve with its 81-bit subgroup.

The exact pair index and batched affine adder passed the n=83 reference-law
tests, including exceptional additions. The probe verified all pair sums
against scalar point additions, found a planted four-summand witness, and
proved the selected outside-range target absent from this small supplied
point list. The point list consists of consecutive generator multiples. It
is an arithmetic control, **not an admissible natural-yield factor base**.

| Points | Unordered pairs | Pair build reference / batch median ms | Pair ratio | No-witness query scalar / batch median ms | Query ratio |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 16 | 136 | 0.544 / 0.181 | 3.00× | 0.572 / 0.078 | 7.36× |
| 32 | 528 | 2.142 / 0.513 | 4.17× | 2.218 / 0.275 | 8.05× |
| 64 | 2,080 | 8.619 / 1.434 | 6.01× | 8.817 / 1.046 | 8.43× |

All three paired timings for each arm are in [raw.jsonl](raw.jsonl).
The implementation was `f6_wide_geometry.rs` SHA-256
`e6c2bf0a2f91a1444efcd1a7c7d1010c6edf27d6493d88dc3f6b93f1d63efb32`;
the probe was `f6_wide_n83_probe.rs` SHA-256
`5a6a604616639ddfbce312bdeca53566bbedb0455c1bc0d54a7827ea10274795`;
the protocol was SHA-256
`d3791cdaf253cb0b10c1a48a78bcfcd5e6ef4df9f1b6ce1d257c5c2bc165663b`.

These are exploratory component timings on an unisolated host. They do not
measure a complete ordinary-query PDP, F6 versus F4/F5, factor-base
construction, relation yield, final LA, target recovery, or rho. The current
F6 Boolean backend still cannot encode the four-summand n=83 chain, and the
pair closure is a known meet-in-the-middle component. No n=83 algorithmic or
end-to-end IC speedup is established by this result.

## Exact coverage ceiling for a four-summand tail

For any base of `B` distinct usable subgroup points, at most
`binomial(B+3,4)` unordered four-point multisets exist. Each has only one
group sum. Thus a uniformly drawn subgroup target has four-summand coverage
at most `min(1, binomial(B+3,4)/r)`, with no independence or distribution
assumption. Repeated sums only lower the coverage. On the tested K_1 subgroup,
`B=4,096` gives at most **0.1371%** coverage and requires 8,390,656 pair
entries; reaching a 1% ceiling needs at least `B=6,733` and 22,670,011
pairs. On the K_0 n=83 confidence-gate subgroup (`r =
2,417,851,639,230,796,216,685,689`), the same `B=4,096` ceiling is
`4.86e-12`; even a 1% ceiling needs `B=872,790`, or 380,881,628,445
pair entries. These are ceilings, not predicted yields. The wide pair path
is useful as an exact terminal closure and arithmetic kernel, but a complete
n=83 F6 design needs a higher-arity strategy or a different relation law.
