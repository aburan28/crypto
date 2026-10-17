# Issued checkout and native-control commands

Working directory for all commands was `/Volumes/SSD990-2/crypto/isogeny-degree-challenge-20261015`, merged source `3ce9694072773badcedc815ddbb953339888e433`.

Before adding the fixture-local `.gitattributes`, the actual Git checkout used:

```sh
git -c core.autocrlf=true -c core.eol=crlf checkout-index \
  --prefix=/private/tmp/isou-fixed-modulus-de577e/checkout-before/ -- \
  research/p192_p224_large_degree_isogenies_20261009/kernel-p224-1471/search.json \
  research/p192_p224_large_degree_isogenies_20261009/kernel-p224-1471/p224/ell-1471.json \
  research/p192_p224_large_degree_isogenies_20261009/kernel-p224-1471/p224/ell-1471.stderr.txt
```

After adding that attribute file, the same command used the fresh prefix `/private/tmp/isou-fixed-modulus-de577e/checkout-after/`. These are historical commands in their before/after worktree state; running the first command on the repaired checkout also preserves bytes. `git-attributes.txt` records text=unset for the three exact paths. `verification.json` checks the after-checkout bytes against both working-tree files and committed Git blobs, without exclusions.

For each arm `reference` and `candidate`, the existing frozen native integration-test executable was invoked as:

```sh
ISOGENY_FIXED_MODULUS_EVIDENCE=/Volumes/SSD990-2/crypto/isogeny-fixed-modulus-builds-20261016/windows-fixture-controls-20261016/ARM \
ISOGENY_CONTROL_BINARY=/Volumes/SSD990-2/crypto/isogeny-fixed-modulus-builds-20261016/clean-retake-v1/ARM/install/isogeny \
/usr/bin/time -l FROZEN_TEST_EXECUTABLE \
  --exact accepted_p224_1471_cache_control_and_rehashed_negatives --nocapture
```

Reference test executable: `clean-retake-v1/reference/target/release/deps/fixed_modulus_replay-df8b92c9bfc1708a` under the build root above. Candidate: `clean-retake-v1/candidate/target/release/deps/fixed_modulus_replay-0b3ab357f28881ab`. Each arm's stdout, stderr/OS receipt and exit code, and every actual judge-child command/input/output/trace/resource receipt, are retained in `native-controls.tar.gz`. Reproduction must use a fresh output directory, not overwrite these historical paths. This test selection is a bounded checkout/cache retake; the earlier complete native suites are separately preserved.
