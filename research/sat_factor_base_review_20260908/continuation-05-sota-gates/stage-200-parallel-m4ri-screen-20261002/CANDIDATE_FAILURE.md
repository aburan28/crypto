# Stage 200 first-candidate failure and scratch correction

The first frozen-order serial control completed successfully. The following
parallel-clear candidate terminated with return code 101 after four Rayon
workers reported:

```text
thread '<unnamed>' panicked at src/cryptanalysis/pq_f4_f2.rs:532:35:
RefCell already borrowed
```

The failed process consumed 22.299680 wall-seconds, 161.557050 total
core-seconds, and 3,357,638,656 bytes peak RSS. It produced no solver JSON, so
it is invalid timing and cannot satisfy any correctness or speed gate. Its
receipt and stderr remain immutable and charged.

The cause is the newly introduced nesting: an outer fixed-X1 F4 call held a
mutable borrow of the worker's thread-local full-M4RI scratch while its
parallel row-clear scope let that worker execute another outer F4 call. The
second call attempted to borrow the same thread-local slot.

The correction takes the scratch value out of its thread-local `RefCell`
before entering elimination and leaves a default value in the slot. A
re-entrant outer call can therefore take independent scratch. Normal completion
returns the original allocation to the originating thread's slot. Pivoting,
tables, rows, thresholds, commands, and all scientific parameters are
unchanged.

The differential test is extended to run twelve outer M4RI calls concurrently,
each forcing nested parallel row clearing, and to require exact rows, pivots,
logical XORs, performed XORs, and table-preparation XORs against the serial
result. After that regression passes, build a fresh corrected binary and
restart the complete frozen pair in the original `serial, parallel` order.
The original serial result is retained for accounting but is not paired with
the corrected candidate.
