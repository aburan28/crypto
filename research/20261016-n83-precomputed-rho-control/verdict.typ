// Verdict under the decision rule frozen in PROTOCOL.md (median paired ratio IC/arm, 95% CI).
*Verdict (decision rule of `PROTOCOL.md`).* Neither online "IC beats rho"
reading survives a matched-precomputation rho. For every time-matched and
memory-matched comparison, the 95% interval of the median paired ratio
IC/precomputed-rho lies entirely above 1: precomputed rho is faster online.

#table(
  columns: (auto, auto, auto, auto, auto),
  align: (left, left, right, left, left),
  table.header[*claim*][*matched arm (run, runner)*][*median IC / rho-with-table*][*95% CI*][*rule*],
  [K3300, 113×], [(a′) = IC wall 680 s, 1 thread (run 2, EPYC 7763)], [3.74], [[2.65, 11.05]], [does not survive],
  [K3300, 113×], [(b) = IC max RSS 10.79 GB (run 2, EPYC 7763)], [40.61], [[24.27, 71.72]], [does not survive],
  [K3300], [(b) (run 1, EPYC 9V45)], [42.72], [[25.49, 74.91]], [does not survive],
  [K600, 12.29×], [(a′) = IC wall 772 s (run 2, Xeon 6973P-C)], [81.52], [[46.35, 236.59]], [does not survive],
  [K600, 12.29×], [(a) 650 s ≤ IC wall 785 s (run 1, EPYC 7763)], [214.73], [[47.94, 483.84]], [does not survive],
  [K600, 12.29×], [(b) = IC max RSS 5.75 GB (run 2 / run 1)], [395.53 / 324.88], [[177.9, 1134.5] / [148.1, 936.2]], [does not survive],
)

Run 1's k3300 arm (a) (650 s on an EPYC 9V45 where IC's setup took 498 s) is
not time-matched and is excluded from the verdict (amendment 1); it reads
16.26 [7.23, 22.60]. Positive control: on the same runners the empty-table rho
is slower than IC, 70.9× [19.2, 92.4] (K3300, run 2, EPYC 7763) and
3.62× [0.93, 4.48] (K600, run 2), so the jobs reproduce the direction of the
original readings. The K600 reading is 3.6–5.3× here against 12.29×
originally because this rho is the fold worker; the hosted run's reference was
the slower strong rung 3.

*Operation counts (AGENTS.md §5a).* Wall time is the claims' metric and the
verdict's. On counted group-addition equivalents no ranking is admitted: K600
IC emits probes and root calls, and the IC/rho order flips between 0 and 0.1
GAE per native operation (IC's counted additions are about 5 per target when
its probes are free); K3300 IC (`koblitz_orbit_dlp_fast.rs`) emits probes but
no root-call count, so its count is a lower bound and cannot be ranked. The
bands are in `analysis.json` (`operation_count_band`).
