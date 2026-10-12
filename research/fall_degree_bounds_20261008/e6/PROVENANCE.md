# E6 binary provenance

- **Part A** ran on a copy of `hky_last_fall_degree` (`sha256_partA`), the same build E4
  used.
- **Part C and the start of Part B** (cell 3:18:6, started 00:52Z) ran on `hky_grid`,
  sha256 `124f121ae8e85dfcf4ac10251eeeeb19361879987c6cfa4fe7b50a2c8d00e099`.
- At 2026-10-12T00:56:28Z, `hky_grid` was rebuilt in place to add the Part D `control`
  subcommand. The new sha256 is
  `1f2b375eefb741f43a858555ab9d499a593b8ecfe216e3c0cda19f3cefeb7b7a`.
- Part B cells launched after that time used the rebuilt binary. The `newcells` and `sat`
  code paths are unchanged in the diff, which adds only a match arm and two imports, so
  their outputs do not depend on which build ran them.
- Part D used the rebuilt binary.
