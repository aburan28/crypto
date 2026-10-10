# Fresh-public-target confirmation: registered before execution

The first panel in `../RESULT.md` found a 181–224× complete-solve F5 versus
inherited-F4 difference, but the `/usr/bin/time -l` wrapper failed after
each worker wrote a complete report. This separate confirmation uses plain
`gtimeout` and records its exit status, with no post-run `sysctl` call.

The new public point is fixed by the worker's hash-to-curve fixture seed
`20261003018`; its scalar is never constructed or supplied. Generate it once
before measurement and freeze the exact coordinates. Use algorithm seed
`20261003033`. The same binary SHA-256
`0b4a97c499ae4b110e4606ae746eb7cd2e47aebd061e8397593c5c5982b02c3a`
and mathematical state SHA-256
`edbff76da6442b9f2e5e8235682c9bf1052465f310c765ba8a37d018a60bf107`
apply. The existing F5 and inherited-F4 candidate IDs remain unchanged:
only the frozen one-target workload and run IDs change. This corrects the
first protocol's suggestion to give an unchanged method a new candidate ID.

One worker, n17a1, 62 usable factor-base points, 29 folded columns, m3,
degree-3 Macaulay limit, node budget 8192, batch trials 1, target cap 32,
dense final LA, imported certified logs, fixed symbolic template and five
exclusive online IC phases remain identical to the first panel. Rho uses
signed Frobenius, one worker, 65,536 maximum iterations per restart and the
same public point. Setup, fixture construction and imported log verification
stay outside the online intervals; each online interval includes scalar
replay. All target-dependent failed attempts stay inside target PDP time.

Run exactly three fresh processes per arm in this fixed order: inherited F4,
F5, rho; F5, inherited F4, rho; inherited F4, F5, rho. Use the same 600-second
cap for every process. Retain stdout, stderr, start/end times, exit status,
source/binary identity, phase ledgers, attempts and all failures. The
predeclared confirmation gate is a median of paired F5/inherited-F4 verified
online ratios at least 2.0×, all nine processes exiting zero, all nine
recovering the same scalar with replay, and all five exclusive IC phases
present and summing to the reported interval. Report the median and range
of paired rho/IC ratios whether they favor IC or rho. No target reselection,
run replacement, or security-scale extrapolation is allowed.
