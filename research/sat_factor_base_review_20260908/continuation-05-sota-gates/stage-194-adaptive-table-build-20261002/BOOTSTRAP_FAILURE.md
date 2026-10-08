# Stage 194 screen bootstrap failure

The first frozen control launch stopped before spawning the backend because the
native meter requires its output directory's parent to exist and
`development/screen` had not yet been created:

```text
refusing to reuse output directory .../development/screen/01-current-r1:
No such file or directory (os error 2)
```

No solver process ran, no timing was produced, and no output directory or
receipt was created. The additive
`development/screen-layout/receipt.json` charges creation of the missing
parent. The unchanged frozen order then ran once: current followed by
adaptive-table-build.

