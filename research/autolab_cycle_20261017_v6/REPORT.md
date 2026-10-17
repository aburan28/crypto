Full binary norm predicate checkpoint

Date: 2026-10-17T03:52:31.358396+00:00

The actual GF(2^83) fixture now compiles the full norm predicate over174 primary bits. This includes source coefficient construction, the exact h section/trace condition, unique k/s recovery and both remaining field residuals. It is not a proof that every original S3 root belongs to this norm chart.

Schoolbook:33524 reachable AND,56594 XOR,327115 standard clauses. Cutoff8:13247 AND,44885 XOR,219448 clauses. AND reduction=60.485%; clause reduction=32.914%. These are representation counts on the same full predicate, not wall-speed or all-roots F4/F5/F6 gains.

Independent literal-norm evaluation passed73216 complete primary assignments across8 small curve/target models, including even extension degree, a0/a1 and c=0. The native search bridge separately recovered all28 primary models of its tiny case. Its actual hard search has NOT been registered or executed.

The backend is a pinned modified CryptoMiniSat functional-primary Gaussian fork, not F4/F5/F6 or unchanged upstream CMS. Native one-thread/configuration readback and model replay are implemented. Library-dependent test is normally ignored in CI but was explicitly run locally. Native library is an external dependency, not uploaded in the source-only PR.

Potential admission issue: total XOR44885 and variables58307 exceed the old per-matrix caps16384/32768. This is not proof that an individual connected matrix exceeds them; inspect components/counters before deciding whether to change the backend. Any UNKNOWN/refusal remains inconclusive and cannot imply an empty root set.

Next: freeze a bounded native hard norm search with full dependency identity, then restore every original ordering/auxiliary fiber and ghost/exception coverage before all-original-root claims. Goal8x, novelty, requested whole field/curve coverage and matched baselines remain open.

Reports remain local because automatic approval review rejected publication without explicit report-payload authorization.
