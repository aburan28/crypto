# Evidence delivery validation

Local validation on macOS ARM64, Python 3.12:

- The final portable replay passes every original SAT, generic and local-native
  audit, verifies source/binary/registration bindings and reconstructs the F5
  nullspace. It never executes an archived producer. `REPLAY-v2.json` retains
  the new audit duration separately from original measurements.
- Four evidence tests cover complete transported replay, original table keys,
  observations and failed rows, a changed archive, and equivalent versus
  genuinely different column-point encodings.
- The complete tournament suite passes **284 tests: 281 pass and three justified
  platform skips**. The process/RSS accounting control runs with host process
  access. Exact command:

  ```sh
  PYTHONPATH=research/ic_candidate_tournament_20260915 python3.12 -m unittest \
    discover -s research/ic_candidate_tournament_20260915 -p 'test_*.py'
  ```

- `git diff --check` passes. The archive is below GitHub's per-file size limit;
  its committed inventory verifies 28,948,695 compressed bytes and 3,851 files.

These checks validate retained evidence and reporting. They do not establish
isolated timing, complete SAT preexecution coverage, a successful F5 solve,
family qualification, a matched speedup or the persistent goal's completion.
Hosted checks and review remain the exact-head merge gates for this PR.
