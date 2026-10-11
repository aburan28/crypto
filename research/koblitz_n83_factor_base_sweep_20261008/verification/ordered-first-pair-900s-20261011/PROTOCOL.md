# Ordered first-pair N83 bounded comparison

The owner authorized at most 900 additional local experiment seconds after
9,900 charged seconds. This gate charges the full declared cap for each
launched stage, including an early exit or infrastructure failure. The staged
caps are 300 seconds for one clean-source Linux worker build, then 300 seconds
each for ordered and unrestricted first-pair search on the same retained
K=2,048 public-x-hash, signed-Frobenius base. The cumulative charge cannot
exceed 10,800 seconds. Tests and code review are preparation; no search starts
until they pass and the source is committed.

Freeze the exact source commit, ELF SHA-256, pinned Linux image ID, replayed
panel object hash, fixture hash, command and cgroup state. The build uses an
offline, network-disabled Linux arm64 container, four CPUs, 8 GiB memory and
a 300-second external wall cap. Search uses the existing
`chain_search_supervisor.py` with one CPU, 4 GiB memory, zero swap, no network,
and a 300-second external wall cap for each run. The supervisor checks input
and source snapshots, worker receipt contracts and container cleanup.

Both search arms use curve a=0, public fixture 0, K=2,048, policy
`public_x_hash`, seed `2026100801`, five summands, one trial, 150,000 SAT
variables, 1,000,000 orbit-domain clauses, 32 models and 1,000 SAT conflicts.
The ordered arm additionally constrains `x0 <= x1`; the second arm has the
same orbit domain without that constraint. Run ordered first. Launch the
second arm only if the first has a valid supervised receipt, cleanup is
verified and its declared 300-second cap remains available. Do not retry a
failed or timed-out launch within this authorization.

This comparison can establish construction counts and a censored search
outcome at the same conflict limit. It cannot determine total cold
index-calculus runtime or a best base without natural relation yield, rank
growth, column and individual logs, and a matched rho run. L0 container wall
times on this shared host are descriptive and do not substantiate a speedup.
