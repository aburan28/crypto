# Retained pre-outcome CI failure

The first independent-verifier freeze at `aab6c4b` had source/hash
preflight PASS but [the archive-replay CI job](https://github.com/aburan28/crypto/actions/runs/36162068734)
failed on a fresh checkout before any selected toy outcome: Git does not
materialize the empty `evidence/` directory, and `ci_replay.py` called
`iterdir()` on it. The corrected pre-outcome head accepts an absent
empty directory as `OUTCOME_HELD`; a partial directory without a receipt
still fails. No toy enumeration, width decision run, or receipt was produced
under the failing head.
