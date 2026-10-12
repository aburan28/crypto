# Bounded correction to the second implementation

The [second comparison](records/attempt-2/comparison.json) passes both
instruction gates (F4 5.5122%, F5 11.0349%) but **fails memory acceptance**.
On the high-bit F4 holdout, additional peak live bytes increase from 54,104
to 54,112, consistently across primary and paired runs. Keep the comparison
and its failed acceptance flags unchanged.

The [native layout probe](records/layout_probe.rs), with its
[recorded output](records/layout-probe.json), measures the original two-way
index enum at 48 bytes, the three-way implementation at 56 bytes, and a
nested fast-index implementation at 48 bytes under the recorded Rust compiler.
The column layout is heap-resident in an `Rc`; the added enum metadata
explains the eight-byte fallback regression. This is a representation defect
in the same memory-bound hypothesis, not a new workload or optimization target.

Correct the representation: preserve the original outer standard/fast enum
and nest hashed/dense storage inside the fast arm. Keep table selection,
the byte budget, column values and all fallback/control semantics unchanged.
Add a unit assertion against the original enum's metadata size, alongside
the independent lookup checks. No extra boxing or allocation is needed.

This explicitly extends AMENDMENT-2's two-attempt stopping rule by **one
bounded correction and complete rerun**. No acceptance threshold or input is
relaxed. Save the corrected results at `records/attempt-2-fixed/`, reusing the
same two frozen baseline profiles and running fresh candidate primary,
five A/B/A′ pairs and the two-flag compatibility control. Retain superseded
validation status accurately. Run the required complete checks on the corrected
source. No further hypothesis or adaptive tuning follows this correction;
if it still fails, revert the direct-index change and deliver the first
candidate's independently accepted F5 improvement with F4 allocation-only
results explicitly labeled instruction-neutral.

Include both rejected attempts and this correction in the final source-linked
report and PDF. Preserve all raw evidence; do not hide the eight-byte regression
inside an aggregate memory value or infer runtime performance from instructions.
