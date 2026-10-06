# Amendment 02: full-universe byte arithmetic

Date frozen: 2026-10-06

This amendment is frozen before source implementation or execution.

Amendment 01 states the correct universe-bounded formula but gives an
incorrect decimal total.  For `U = 34,562,148,612`:

```text
ceil(U / 64) * 8 =   4,320,268,584 bitmap bytes
16U                = 552,994,377,792 value bytes
stream accumulator =              16 bytes
total              = 557,314,646,392 bytes
```

The corrected total is about 557.31 GB (decimal), or `2^39.019` bytes.  The
earlier 557,480,271,372-byte statement is withdrawn.  No experimental result
was produced under the incorrect decimal total.
