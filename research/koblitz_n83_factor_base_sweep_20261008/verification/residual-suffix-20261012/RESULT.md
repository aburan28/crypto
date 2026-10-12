# Checked residual-suffix source gate

The large-prime module now constructs a residual signed-Frobenius base from
the exact ordered suffix of two replayed, prefix-related bases. The degree-7
control extracts the intended suffix and rejects changed prefix points and
changed suffix orbit labels. This is a source and small-curve correctness
check; no retained N83 suffix construction, partial search, relation yield,
graph completion or cold-runtime measurement was run.

The focused release module passed 4/4 tests. The full release library suite
passed 3,986 tests with 112 ignored, and the study Python suite passed 39/39.
The touched Rust file passed `rustfmt --check` and the source diff passed
`git diff --check`. The repository-wide format check encounters missing bench
files and unrelated formatting differences in the current checkout; it is
not a validation result for this patch. The [machine receipt](validation.json)
binds the changed Rust file's SHA-256 and three compressed raw test logs.
