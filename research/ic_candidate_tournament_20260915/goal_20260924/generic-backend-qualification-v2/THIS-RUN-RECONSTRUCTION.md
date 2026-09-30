# Independent reconstruction of seed 2026092902 exposures

The censored second registration
([Actions run 36580669479](https://github.com/aburan28/crypto/actions/runs/36580669479),
seed `2026092902`) lost its campaign archive. Treat every public point that
dispatch's prepare path could have generated as exposed, whether or not a
timed-out worker reached it.

[this-run-exposures.json](this-run-exposures.json) freezes that corpus
(SHA-256 `64f30e19f0ef1a4c4b95d4168c977b13de9559cb8a934cfcc9c3bf8a63fa0a25`).
It was produced by
[reconstruct_this_run_exposures.py](reconstruct_this_run_exposures.py) against
the pinned worker `765c3c5f19032bd852163805f257c56babef2040`. The generator
excludes the sealed three-round history, the seven supplemental fixture
corpora, and the first-run
[lost-campaign-exposures.json](lost-campaign-exposures.json)
(SHA-256 `a728677b199eac02800d8338204d5306f391ec5da757c910bec1e51955fe7b41`).
All 25 accepted points are distinct from that first-run set. A second local
reconstruction matched the checked file byte for byte.

This is an exclusion census only. It does not recover measurement receipts,
natural yield, family admission or competitive costs from the censored run.
Never redispatch seed `2026092902`. A later registration must pass both
exposure corpora to `--exposed-fixtures` before sampling fresh targets.
