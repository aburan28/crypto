# Native metadata regeneration

The experiment, analysis, replay supplement and dashboard additions use Rust.
The renderer only draws and lays out recorded numbers. No legacy Python
research generator is extended or run on the new sessions.

ecbench_open_questions_publish rebuilds the report, per-session Markdown,
render-data, requirement matrix, scoreboard panel, leaderboard HTML/JSON,
diagnostic timeline series, and lab-browser curve/method/factor-base/session/
yield indexes. It adds immutable records and native summaries to the exact
frozen metadata at f9011e5375155bd9ef7839abe4eea12718e4b401. Unchanged
historical JSON spans retain their bytes; historical root arrays are preserved.
No historical mathematical analysis is recomputed. Count quotients remain
diagnostics.

The native-publish crate compiles that same source for light metadata CI with
standard SHA-256. Its --check compares generated outputs byte-for-byte.
Historical generator checks remain in a temporary baseline worktree; a source
hash guard checks that their inputs have not changed. New outputs use the
native checker. Those historical checks run no new research. Later historical
input changes require a native regeneration/version update.

Run from the repo root:

    cargo run --release --locked --manifest-path research/ecbench_open_questions_20261014/native-publish/Cargo.toml -- --check

Never overwrite a frozen session. Later sessions require new directories.
