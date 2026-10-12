# Stage 199 unset-replay bootstrap failure

The first unset replay stopped before spawning the backend because the charged
layout command had already created `development/default-replay`, while the
native meter requires its output directory not to exist:

```text
refusing to reuse output directory .../development/default-replay:
File exists (os error 17)
```

No solver process ran and no timing or receipt was produced. The corrected
immutable layout uses `development/default-replay/run` as the meter receipt
directory and stores composition outputs inside that run directory. The unset
backend command and every scientific parameter remain unchanged.

