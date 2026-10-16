# Preserve exact public fixture bytes on Windows

The failed Windows CI job and actual Git autocrlf reproductions are retained here. See verification.json for the exact paths, SHA-256 values, native test arms and zero-exclusion output comparison. native-controls.tar.gz contains all 106 raw files; native-manifest.json binds every member. ../custody-approved-20261016/README.md and report.pdf give the combined source/custody context.

The fixture-local .gitattributes changes checkout treatment only. No fixture bytes, expected hash, arithmetic, assertion, or admission requirement changed. Windows CI on the repaired commit remains pending.

Exact issued commands and before/after checkout state are recorded in [commands.md](commands.md).
