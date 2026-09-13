# Stage 108: direct-routed support-table selection

Stage 108 archives the measured Stage 105--107 optimization chain. The n=53
factor base remains the algebraically point-defined 9,964-point, 94-orbit base.
Its construction does not enumerate the target subgroup or use discrete-log
labels. Each complete run derives all 94 base logs and the target log from 95
relations that are checked in the group.

Stage 105 partitions the 24,805,379-entry exact support table into four equal
hash ranges. Expansion and insertion can then run without locks on separate
shards. On the identified Intel Xeon 6973P-C host, the pipelined sharded arm
reduced process wall from 3.994782 to 3.302822 seconds, a ratio of 0.826784.
Support setup improved from 2.272968 to 1.533798 seconds, or 1.481922 times.
The arm increased direct core-seconds by 25.95 percent and did not beat rho on
that run, so sharding alone was retained only as an optimization candidate.

Stage 106 replaces the splitmix shard-selection hash used on every query with
the xor of low and high x-coordinate windows. Five independent, host-identified
runs all favor direct routing over mixed routing. Their wall ratios are
0.951938, 0.899365, 0.870967, 0.904540, and 0.906388, with median 0.904540.
The median direct/mixed core ratio is 0.912199. The same direct arm also beats
the matched automorphism-optimized rho process in all five runs, with ratios
0.975768, 0.860651, 0.915716, 0.831349, and 0.905680. The median is 0.905680.
The receipts cover AMD EPYC 7763, AMD EPYC 9V74, and Intel Xeon Platinum 8573C
hosts and include self-hashed CPU, cache, affinity, kernel, and feature data.

Stage 107 selects the direct-routed four-shard table in the four-core known and
unknown target launchers. The public-hash unknown-scalar run on EPYC 7763
recovers scalar 7,892,094,459,170 without constructing or supplying it to the
collector. Direct wall is 4.307977 seconds and rho wall is 4.362556 seconds,
for a ratio of 0.987489. The recovered scalar is accepted only after its public
group equation reconstructs the target.

The selected five-pair panel ran on an AMD EPYC 9V74 host. Its direct/rho wall
ratios are 0.847095, 0.858937, 0.839190, 0.829081, and 0.838930. Direct wins all
five, the paired median is 0.839190, median direct wall is 3.630413 seconds,
and median rho wall is 4.308264 seconds.

The full-cost boundary remains unfavorable. Fresh build plus median selected
direct is 17.636692 times one rho target in the Stage 107 panel. Direct uses
roughly 1 GB peak RSS and more core-seconds than rho even though its online wall
time is lower. Licensed same-instance Magma F4, the complete growing-n solver
panel, and unaffiliated reproduction and novelty review remain absent. These
receipts establish a repeatable public-synthetic online engineering crossover,
not a full-cost or Koblitz index-calculus SOTA result.
