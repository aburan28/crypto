# Integration checks after publication

The first PR #350 head, `7259f1f6593078bb5d44299b0169b2079b0c0758`,
contains the measured implementation and complete accepted evidence. Its
benchmark source hashes remain the measurement identities. The subsequent
integration fixes do not modify those solver sources or measured inputs.

GitHub checks exposed these compatibility problems:

- The polynomial workflow's unquoted commands ending in `::` were invalid
  YAML. Quote those commands and install the compatible frozen dependency
  lock used by the local measurements.
- Three current-code workflows copied the older Stage-20 lock, which predates
  the optional Redis dependency. Their operational smoke builds now use the
  compatible frozen lock. The Stage-21 and Stage-23 runners still require the
  archived lock for production; the smoke lock is explicitly recorded and
  verified, and smoke evidence remains ineligible for scientific promotion.
  Fetch the separately pinned Stage-24 helper dependencies before its offline
  build, since its archived `cc` version differs from the current root lock.
- Two older examples explicitly initialized every index-calculus option but
  omitted `weil_charts`. Set it to `None` for their existing SAT strategy.
- The copied Pages scoreboard has no adjacent research directory. Link its
  thirteen evidence references to immutable files at the published evidence
  commit. The repository HTML remains the canonical file copied by the site
  builder.

Local checks: all examples and tests type-check with the compatible lock;
all 17 site-build tests, 13 Stage-21 custody tests and 12 Stage-23 custody
tests pass. All eight Rust producer tests pass. New lock tests preserve the production pin and reject unbound
workspace locks. Build/check logs, including the initial failures, are kept
in this directory. GitHub CI is rerun on the integration commit.

Initial failing runs: [site](https://github.com/aburan28/crypto/actions/runs/34942858826),
[Stage 21](https://github.com/aburan28/crypto/actions/runs/34942858566),
[Stage 23](https://github.com/aburan28/crypto/actions/runs/34942858569),
[Stage 38](https://github.com/aburan28/crypto/actions/runs/34942858622),
[SOTA reproduction](https://github.com/aburan28/crypto/actions/runs/34942859229),
[workflow syntax](https://github.com/aburan28/crypto/actions/runs/34942845669).
