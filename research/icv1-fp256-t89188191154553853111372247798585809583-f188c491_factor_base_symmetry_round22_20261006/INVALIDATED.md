# Superseded round-22 artifacts

`symmetry-result-preprojection-fix.json` is retained because research runs are
not overwritten.  Its exact factor-base rebuilds, point manifests, and
transport censuses agree with the canonical artifact, but its
`projected_138031_rows_log2_operations` field multiplied by 138,031 after the
cross-colour birthday formula had already priced the requested collision
count.  Its projection block is invalid.

`symmetry-result-pre-la-fix.json` corrects that double count but predates the
explicit sparse-Wiedemann and Berlekamp--Massey lower-bound field.  Its
projection is incomplete and superseded.

Only `symmetry-result.json`, SHA-256
`3116c678d257794040c5519c85a8037177e12381e5f948a47cf1335d0cafde16`,
is canonical.
