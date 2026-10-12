# Matched retained-base N83 search gate

The owner authorized **at most 900 additional local experiment seconds** after
10,800 conservatively charged seconds. Freeze this protocol before either
search. Charge each launched stage its full declared cap, even on early exit,
timeout, or producer failure: 280 seconds for one source-matched offline Linux
build, 280 seconds for each of two searches, and 60 seconds reserved for
cleanup and receipts. Do not retry a failed stage within this grant. Total
conservative charge is at most 11,700 seconds.

Compare the retained primary `a=0`, `public_x_hash`, seed `2026100801`,
K=1,182 and K=2,048 signed-Frobenius bases. They are exact prefix-related
objects, verified byte for byte in `verification/v2-prefix-lp-20261012/`.
Use the same public fixture-zero target, five summands, ordinary unrestricted
Frobenius-orbit-phase domain, one trial, at most 150,000 variables, 1,000,000
domain clauses, 32 models, and 5,000 SAT conflicts **for each** search.
Run K=1,182 first; run K=2,048 only after the first has a valid supervisor
receipt and verified container cleanup. A solver cap is censored, not a
relation-yield estimate.

The existing `chain_search_supervisor.py` must freeze exact source commit and
hashes, Linux ELF hash, replayed object and fixture hashes, commands, cgroup
state, output, and cleanup. Use its immutable local Linux arm64 image ID,
one CPU, 4 GiB hard memory, zero swap, no network, and 280-second outer wall
cap for each search. Build a worker from the clean committed source with the
existing `ordered-first-pair-900s-20261011/run_build.py` only as a build
wrapper: its receipt schema names its original study, but the resulting
binary has no ordered-pair setting until the search command requests one.
Build under its pinned Linux image, four CPUs, 8 GiB hard memory, no network,
and 280-second outer wall cap. Reuse neither the old binary nor an unverified
object. Locally retained compressed objects must match the frozen SHA-256 and
committed manifest/replay/upload receipts before launching.

Record both exact search commands and full supervisor output. Independently
audit source, binary, object, target, cgroup, counters, and cleanup. Report
every setup, model, relation, rank, linear-algebra, and individual-log stage
separately. These one-trial capped searches cannot select a minimum complete
index-calculus runtime factor base, even if their outer wall times differ.
