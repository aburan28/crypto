# Charged rho retry amendment, before batch rerun

The frozen V1 command for `(p,r,round,index)=(61,19,0,0)` exited 101. It
completed field and pairing-input selection and wrote one setup record to
[`RUN_61_R0.jsonl`](RUN_61_R0.jsonl). Its first same-target rho call did not
return a verified scalar under seed `0x4241544348001000`; the source asserted
and stopped before a target record. The [stderr](RUN_61_R0.stderr) retains the
exact panic. This is a reference-walk retry condition, not a mathematical
result or a zero-yield pairing row. Do not overwrite this run.

For V2, keep the public target batches, H search, field, pairing algorithm,
finite-field table, three rounds, and per-attempt rho cap unchanged. On each
target, attempt rho seeds
`seed_base + retry * 0x01000000 (mod 2^64)`, for `retry=0,...,15`, where
`seed_base` is the V1 target seed. Record each attempt's seed, result, CPU,
wall, and field counts, and charge every failed attempt in the cumulative rho
reference. Stop with an explicit censored row if all 16 attempts fail; do not
substitute a successful later seed without counting earlier attempts. A
pairing/rho ratio is admissible only for a prefix in which every rho target
completed and every pairing scalar independently replayed.

This amendment is versioned before the V2 run. It changes the rho stopping
rule only; no target, pairing input, answer, or selection receipt changed.
