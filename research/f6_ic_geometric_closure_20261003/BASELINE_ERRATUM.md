# Earlier inherited-F4 manifest split-rule field

The immutable `research/ic_solver_online_20261003/inherited_f4-record.json`
records `point_decomposition.limits.split_rule: "lowest-free"`. The exact
worker and Gröbner source used by that panel instead call
`split_rule_default()` with no `SOLVER_SPLIT_RULE` override, resolve its
`Auto` value under `SolverEngine::InheritedF4`, and select `HighestFree`.
`MatrixF5` resolves `Auto` to `LowestFree`.

The earlier panel's raw points, verification, intervals and solver-family
comparison remain retained. The inherited-F4 candidate identity hashed an
inaccurate descriptive field, so that ID must not be reused as the identity
of a correctly specified baseline in this F6-IC comparison. The new
baseline record will say `"highest-free"`, carry the new source snapshot and
receive a new candidate ID. This erratum does not retroactively rewrite the
frozen record or promote its Mac timing through the newer CPU isolation gate.
