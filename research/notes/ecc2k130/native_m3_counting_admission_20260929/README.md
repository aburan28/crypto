# Replay the native m3 counting admission

The versioned protocol, source/input lock, exact result, independent replay
receipt, and scoped decision are in this directory. It is a deterministic
integer bound with no random seed or n=131 solver invocation.

From the repository root, replay without overwriting the archived files:

```sh
python3 research/notes/ecc2k130/native_m3_counting_admission_20260929/compute.py \
  --out /tmp/native-m3-counting-result.json
cmp /tmp/native-m3-counting-result.json \
  research/notes/ecc2k130/native_m3_counting_admission_20260929/RESULT.json
python3 research/notes/ecc2k130/native_m3_counting_admission_20260929/verify.py \
  --result /tmp/native-m3-counting-result.json \
  --out /tmp/native-m3-counting-verify.json
```

`VERIFY.json` records the original verifier invocation, including the
machine-specific Python and parent-admission paths. Its semantic fields
(classification, result hash, checks and parent admission stdout hash) are
portable and compared in CI. The committed `RESULT.json` is byte-reproducible.
