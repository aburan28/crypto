# Attribute the inherited F6 basis-build residual

Registered before implementation or timing. The support-local profile
in PR #1465 found that generating rows, collecting columns and packing
account for only 4.75–4.98% of the F4 build timer on prepared n17 T7.
The remaining time includes root basis reduction, child-system
preparation, and inherited-basis specialisation. This opt-in diagnostic
times those three paths, including counts, as nested components of the
existing build timer. It changes no matrix, row space, solver decision,
or default configuration.

Freeze one new exact IC1 candidate on the same prepared n17 Koblitz
curve, 62 actual usable base points, 29 folded columns, archived public
targets T1/T7, imported certified logs, three summands, degree three,
8,192 node budget, 32 target trials, default algorithm environment and
one Rayon thread. Keep all prior optimization flags disabled. Enable
both the support-local and inherited-basis timers, charge their cost
inside the target online interval. Commit source, binary, candidate
manifest, inputs and hashes before timing. Run two fresh processes per
target, preserving every status, stdout/stderr, timestamps, five
exclusive online phases, replay certificate, attempts, reductions,
matrix dimensions and word operations.

Before timing, test that disabling/enabling inherited-basis profiling
preserves an exact prepared F6 target result and that the default-false
flag is omitted from generic `effective_config`. Run the existing F6
geometric-closure control. Require all successful rows to replay the
archived scalar and sum five phases exactly to online wall. If any run
fails, retain it and do not call the comparison a success. Record root
count/time, child preparation count/time, specialisation count/time and
the earlier support-local component times. Compute the unassigned F4
build residual, without double counting nested support-local time.

Decision rule: if specialisation takes at least half of F4 build in
both T7 repetitions, inspect its row rewrite/reduction path next. If
root construction takes at least half in both, optimize root echelon
reduction. Otherwise investigate the largest measured component and
residual before another optimization. Report observed ranges from the
two repetitions and all raw failures. This is a diagnostic, not a
speedup claim. The unisolated Mac cannot promote CPU ratios; n83
ordinary relation yield and one-target IC-versus-rho remain outside
this experiment.
