# Marked CMS build attempts

The first guarded build attempt after the full-rank SAT preparation audit
created `/Volumes/SSD990/llm/tmp/ic-marked-cms-build-20261006-v2` and stopped
in its `extract` phase with `BUILD_FAILED`. The original
[terminal receipt](failed-build-v1/terminal.json) is retained in this repo.
Its source archive matched the
declared SHA-256, but the script rejected `cms/cadical.tar` before applying the
patch, invoking CMake or launching a solver. The failed output and terminal
receipt remain unchanged at that path.

The script and new Rust transport auditor had transcribed the CaDiCaL digest
as `8264713f3dc1c4455162d2912238712bd8030fceab0f4b430d106b5a58058d`.
The extracted bytes hash to
`8264713f3dc1c4455162d2912238712bd8030fceabec0f4b430d106b5a58058d`.
That latter value is also the exact `cms/cadical.tar` digest in the accepted
`static-sat-runtime-v3/native-inputs-macos-arm64/manifest.json` and the native
assets receipt in the separately frozen ordinary registration. Only the two
transcribed constants were corrected. Any later attempt uses a new output
directory; the first failure is not reclassified as a successful build.
