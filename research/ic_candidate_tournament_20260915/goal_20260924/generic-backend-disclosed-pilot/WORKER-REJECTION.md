# Schema-v2 pilot stopped before algebraic work

The [complete compact raw record](RESULT-v2.json) shows that both stage-A
worker processes exited with code 2. Their exact reports say
`undeclared algorithm environment override: KIC_F5_AVX512_UNPACK`.
The runner had set that variable to `0`, while the pinned exclusive-phase
worker rejects **any** `KIC_*` environment override. Each process ended in
about 120 ms with no factor-base construction, natural query, PDP call,
relation, or target solve. The predeclared gate kept all eight stage-B jobs
unexecuted. This is a runner configuration failure, not evidence that F4 or
F5 works or fails.

The schema-v2 jobs are not retried. Before further execution, schema v3 uses
the five *different*, previously disclosed smoke-stage public points and a
new algorithm seed. It removes the undeclared override by giving the worker
only ordinary process variables and `RAYON_NUM_THREADS=1`; it retains the
same source, algebraic configuration, 8 GiB RSS threshold, one-query limit,
60-second cap, stage gate, and independent audit. The prior schema-v1 and v2
registrations and outcomes remain in Git history and in the linked raw record.
