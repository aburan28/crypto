# Stage 193 verifier path correction

The first composed result passed its absolute-path replay `21/21` with
SHA-256
`a7748cb78f5bea886d3374a6d6ad90d4af8ac48c7392f3e4f293e52e8346365c`.
A deliberate repository-relative invocation then failed before producing an
output with:

```text
WDSat config-install command differs from protocol
```

The verifier found the stage root from the relative result path but compared
that relative root against absolute paths in the authenticated meter command.
The experiment, solver receipts, hashes, and terminal interpretations were
unchanged. The correction canonicalizes the discovered stage root before every
artifact and command replay. The first result and verification are preserved
under `development/superseded-relative-path-verifier/`; the preservation
commands and corrected build/test/replay are charged additively.
