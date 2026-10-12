# Stage 103: identified-host n=53 result

Stage 103 archives and replays the hosted Stage 100--102 chain. The factor base
remains the algebraically defined 9,964-point, 94-orbit base. Its selection uses
no target-subgroup enumeration or discrete-log labels, and every run derives
the 94 base logs and target log from 95 verified relations.

Stage 100 replaces repeated binary scalar multiplication in exact reference
validation with doubling tables for the fixed bases `G` and `Q`. It retains the
same reference curve arithmetic, relation hashes, solution, target
reconstruction, all 94 factor-base log checks, all linear rows, and all 95
independent relation replays. The candidate reduces median process wall from
5.235722 to 4.790938 seconds, a ratio of 0.915048. Solution validation is
3.051731 times faster and independent relation replay is 2.817622 times faster.

Stage 101 demonstrates that an anonymous `ubuntu-latest` label is insufficient
for a portable rho comparison. Direct has a 4.488816-second median while rho has
a 3.306672-second median, producing zero direct wins and a 1.353855 median
ratio.

Stage 102 binds the full comparison to this self-hashed host identity:

- AMD EPYC 7763 64-Core Processor;
- family 25, model 1, stepping 1, microcode `0xffffffff`;
- four effective logical CPUs in the process affinity;
- 32 KiB L1 data, 512 KiB L2 per shared pair, and 32 MiB shared L3;
- PCLMULQDQ and VPCLMULQDQ exposed;
- Linux 6.17.0-1022-azure;
- identity SHA-256
  `325e41e6d926514b9843ca08bafa2673f9803ac1bf990667be819cbda2613485`.

Five separately metered direct/rho ratios on that host and the same public
synthetic target are 0.903291, 0.913452, 0.898292, 0.899840, and 0.878913.
Direct wins all five. Median direct wall is 4.628149 seconds, median rho wall is
5.130762 seconds, and the paired median is 0.899840. This passes the
predeclared online threshold of 0.95 for this recorded host class.

Direct peak RSS is approximately 1.055--1.058 GB. Its direct-process core cost
remains 1.58--1.68 times rho. Fresh build plus median direct remains 19.149871
times rho. Licensed same-instance Magma F4, a cross-host-class full-cost
crossover, refreshed selected unknown-scalar and single-core measurements, and
unaffiliated reproduction and novelty review remain incomplete. Stage 102 is a
finite identified-host online crossover, not a full-cost or Koblitz
index-calculus SOTA result.
