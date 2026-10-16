# N83 K3300 root11 Bloom prefilter screen

This bounded screen tests the existing optional three-hash root Bloom filter
against the exact five-byte, fixed `2^30`-slot K3300/root11 index. The only
controlled change is `KIC_ROOT_BLOOM=0` versus `1`. Both arms retain exact
five-byte table lookup after a positive Bloom answer; the table keeps the
first inserted exact-root witness. The Bloom filter is a query prefilter,
with a predicted `2^29`-byte allocation. Both arms run in fresh processes,
in off then on order, using four rank threads, one target thread, rank seed
29, a 1,200-second child cap and a 14 GiB sampled RSS guard.

The curve is `EC1N83Ckb1a1haacb69ee9f37`, polynomial-basis `GF(2^83)`
with modulus bitmask `9671406556917033397649557`, model
`y²+xy=x³+x²+1`, subgroup order `8569786107849059`, and cofactor
`1128547018`. The factor base is the frozen first 3,300 signed-Frobenius
subgroup orbits: 547,800 distinct usable points, 3,300 folded columns,
base hash `a4ba33098d2776520b98a326585a2301f8ccd5bc005df7c1a70fbe8496c61d8b`,
and 903,852,075 exact root keys. A different base, table capacity, layout,
rank or key count fails the comparison gate.

Both arms use the same 1,000 public points from
`research/n83_factor_base_refine_20261014/targets.jsonl`, SHA-256
`35e34151dcc1ad6048939560adf28fc7c1b1017833fe005a656c3f4e2bc07de8`.
The scalar answer sidecar SHA-256 is
`f8cbc90f7f592aae5014ca4f2693a3327a72e935970f22e8f02868e4aeed7f71`.
These points were used in an earlier factor-base panel, so this is a
secondary exploratory algorithm screen. A positive arm needs fresh
same-point Tier 1 candidate/rho pairs for a primary online result.

Compute the paired difference `bloom_on_ms(i) - bloom_off_ms(i)` for all
1,000 points, partition them into 20 contiguous blocks of 50, resample
those blocks 10,000 times with seed `2026101507`, and take sorted 2.5%
and 97.5% endpoints. The Bloom arm advances only if both arms pass full
base, four-point relation and scalar replay, remain within the RSS guard,
and the upper interval endpoint is below zero. Otherwise the selected
K3300/root11 index remains the control. Timeouts, memory caps and verifier
failures remain recorded statuses; they do not become a zero-cost arm.

Record source, binary, fixture and host hashes, complete setup and process
costs, peak RSS, table/filter bytes, root-key count, every target online
phase and probe count. The online interval starts after base/index/log
preparation and includes target query, point decomposition, relation check,
descent and recovery check. The independent verifier replays the complete
factor base and every target scalar and four-point relation.
