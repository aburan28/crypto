# Static zero blocks in the existing polynomial-square table

The sigma-only square-table panel is terminally negative at paired geometric
mean0.996538. Its shared implementation performs65 word loads and65 XORs per
lambda square. This native-only study asks whether any whole `(window,word)`
block is identically zero across all32 entries and can therefore be removed
without a data-dependent branch or a changed field representation.

Freeze the existing `fillSquareTable131`/`squarePolynomialTable131` layout:
13 windows,32 entries,5 output words,8320 bytes. Reconstruct the same table
from the current tracked native arithmetic; list every block's32-entry OR.
Build one static active-block mask from the exact table bytes. The candidate
oracle skips only all-zero blocks, retaining the same low-coefficient spread,
entry indexing, nonzero blocks, and table footprint.

Reference: the current65-load table square. Exact control cases are all131
basis coefficients, zero/one/all-ones, and20,000 deterministic dense canonical
field elements. The masked square must equal both the existing table square
and the independent scalar field reference on every case. All65 blocks must
also reproduce from the independent basis-square construction. No GPU or
whole-walk measurement is authorized by this static protocol.

Admission: at least one whole zero block, exact control PASS, and an exact
positive reduction in loads/XORs. If no block is zero, stop with a negative
static result and do not implement a GPU variant. Even a passing static screen
does not establish throughput: native lowering, resources and a separately
frozen same-GPU whole-walk comparison would be required. Dynamic skipping of
entry0, new packing, ALU-square changes, table shrinking and other map changes
are outside this study. The26 B/s objective remains unachieved.
