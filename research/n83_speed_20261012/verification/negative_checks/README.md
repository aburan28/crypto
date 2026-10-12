# Controlled rejection checks

These are intentional validation tests, not campaign processes. All four
rejection categories returned exit 1. The genuine 139-process panel and frozen
fixtures were never modified.

| Test | Preserved observed rejection |
| --- | --- |
| Copy fresh fixtures; replace training_1 bytes with training_2 bytes | `Fresh point digest mismatch` |
| Pass an all-zero expected SHA-256 for the frozen protocol | `Frozen protocol digest mismatch` |
| Generate into the existing frozen fixture directory | `Refusing existing output directory` |
| Change one copied producer's recovered scalar to `1` | First: `Execution producer digest mismatch`; after changing only the copied execution's producer hash: `Independent native scalar replay failed` |

The copied producer and execution are `wrong_scalar_producer.json` and
`wrong_scalar_execution.json`; they are synthetic mutations of AA01-a used to
exercise failure gates. No solver was launched for these tests. Their original
argv references are archival runtime paths. The negative checker invocation used
a separate one-row panel with expected count 1, outside `/workspace/n83_speed/panel`.
Both wrong-scalar failure receipts are retained. A read-only cloned execution
file initially refused copying; only that test clone was made writable before
updating its digest. No genuine execution receipt was made writable.

Native tools and flags match the adjacent verification README. The deliberate
point-file mutation uses a copied fixture directory; the generator's existing-path
test relies on exclusive creation. These tests support fail-closed validation of
wrong inputs, byte bindings and answers; they do not supply performance evidence.
