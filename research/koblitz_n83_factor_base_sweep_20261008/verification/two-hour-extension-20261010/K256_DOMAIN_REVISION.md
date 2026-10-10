# K256 domain-clause admission revision

Frozen after the first K256 capacity receipt, before any revised launch.
The 3,000,000-clause gate returned `UNKNOWN_domain_clause_cap`; its exact
domain preflight requires 7,195,450 clauses for 21,248 legal x coordinates.
The original outcome is retained and is not a solver verdict. This new fixed
resource hypothesis admits at most 8,000,000 finite-domain clauses while
retaining the same 150,000-variable limit, one CPU, 4 GiB zero-swap cgroup,
no network, K256 stored object, public fixture, source revision, one-model
and 10,000-conflict limits. The extra clause allowance is motivated by the
measured preflight count and was not changed during the earlier run.

Run one new construction gate under 60 seconds. Only on construction PASS,
run a separate zero-trial preflight under 60 seconds. Only if both pass and
at least 620 seconds remain under the two-hour-total allowance, run one
trial with 590-second worker and 600-second outer caps. If resource or wall
limits intervene, retain UNKNOWN and stop. The clean Docker query criterion
from K256_GUARD_REVISION.md remains mandatory. No runtime ranking follows.
