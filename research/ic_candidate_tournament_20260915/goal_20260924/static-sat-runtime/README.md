# Historical static SAT runtime replay

This source bundle replays the v1 and v2 static SAT registrations with their
historical code after the tournament selector changes. It does not rewrite
their source manifests, candidate IDs, workload IDs, panels, or seals.

`freeze_static_sat_runtime.py` recovered each registered Python file from the
Git commit recorded in that run's `source-binding.json` and required its bytes
to match the registered SHA-256. The deterministic archive and receipt retain
36 source files. This recovery happened after execution; it is not a claim
that this archive existed before either measurement.

The original source walker omitted two imported package files:
`producer/evidence.py` and `producer/timing.py`. Both are recovered from the
recorded commits and explicitly listed as **supplemental** files, outside the
unchanged registered manifest. Thus the original manifests do not establish
complete preexecution coverage of the Python runtime. The replay reports
`complete_preexecution_python_manifest: false` and
`promotion_eligible: false`. Passing the replay establishes historical
identity and input reproducibility, not a fully source-bound performance
claim. The scalar and relation correctness audits are separate evidence.

`frozen_sat_runtime.replay()` checks the archive and all source hashes before
executing registration replay in a temporary directory. All local Python
imports come from the archive. Only the existing registered data directories
are linked into that directory. The replay does not launch SAT solvers,
execute a measured workload, or supply live Python code as a fallback.

Live v1/v2 runners still require their exact registered manifest. Tests retain
that failure gate. Future executions using changed runtime sources or a
complete package-aware source walker require a new registration and candidate
identity. The omitted files cannot be retroactively inserted into an old
candidate's identity.

Run historical replay and its failure-gate controls with:

```sh
python3.12 -m unittest discover -s research/ic_candidate_tournament_20260915 -p 'test_static_sat_full*registration.py' -v
python3.12 -m unittest discover -s research/ic_candidate_tournament_20260915 -p 'test_frozen_sat_runtime.py' -v
```
