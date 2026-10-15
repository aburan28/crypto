# V2 shell-wrapper failure before the batch comparison

The first invocation of the committed `run_v2.zsh` exited **1** after
the first source arm completed. In zsh, `status` is a read-only shell
parameter. The script attempted to assign `status=0` before writing
that arm's `.exit` receipt or appending its stdout to `RUN_V2.jsonl`.
The exact [wrapper stderr](RUN_SCRIPT_V2.stderr) says
`read-only variable: status`; its SHA-256 is
`decace5490f6df7e3aa8a409d7c6c3bd4f29e675bd9283d830b7e61791f7f40d`.
The merged [RUN_V2.jsonl](RUN_V2.jsonl) is zero bytes. The first arm's
complete [raw stdout](raw-v2/round-0-source.stdout), SHA-256
`46d48fc6bb6da5acbd7b86599eda03f18b265fac4652522817d23af73d345d39`,
is retained with its empty [stderr](raw-v2/round-0-source.stderr).
That raw record reports 60 verified target rows, but it has not passed
the separate 21-arm batch replay and is **not** used in a comparison.
The wrapper exit happened before it wrote a per-arm exit file; do not
infer one from the raw result. No other V2 arm ran.

The [V3 script](run_v3.zsh) changes the shell variable to
`run_exit_code`, uses new `raw-v3/` and `RUN_V3.jsonl` paths, and refuses
to overwrite them on retry. [RUN_SPEC_V3.json](RUN_SPEC_V3.json) pins
that script's hash before the V3 run. Source, binary and public/checker
inputs are byte-identical to V2. This was an orchestration failure,
not a relation or DLP result.
