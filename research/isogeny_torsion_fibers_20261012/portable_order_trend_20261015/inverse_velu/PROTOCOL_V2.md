# Direct inverse-Vélu fiber oracle, split schedule correction V2

V1 remains frozen as [`PROTOCOL_V1.md`](PROTOCOL_V1.md). Its first
`p=1237,r=1201`, round-zero, one-target smoke run exited 101 in the
inverse polynomial split before producing a child record. The exact
V1 binary hash is in `BUILD_SHA256_V1.txt`; the first stderr and the
independently captured retry exit/stderr are retained as
`SMOKE_R0_N1*`. The experiment did not reach relation collection.

The V1 shifts `c=0..31` were implemented as `[c,0]` in `F_(p²)`.
For three roots in `F_p`, every nonzero `root+c` lies in `F_p` and is
a square in `F_(p²)`. The quadratic-character partition can therefore
be constant for every base-field shift; extending the same V1 shifts
would not discriminate that fiber. This is an algorithmic split-schedule
defect, not evidence against the inverse polynomial identity or the
existence of a fiber.

**V2 change, frozen before its run:** keep all V1 inputs, controls,
caps, 60-child panel, accounting and acceptance gates. Replace the
32 split shifts with exactly `[i,1]`, `i=0..31`, in the fixed field
basis. A failed split still stops the child. Count each attempt and
preserve any new failure. Do not silently broaden to another shift
family or source enumeration. The output schema is tagged V2, with
the source and binary hashes recorded before execution. V1 results
remain a failed precursor and must not be overwritten.
