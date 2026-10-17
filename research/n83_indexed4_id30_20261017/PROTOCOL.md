# Four-byte ID30 exact root index for larger N83 bases

The opt-in `KIC_INDEXED4_ID30=1` with `KIC_INDEXED4_LOAD75=1` stores an
exact-state reference in each four-byte root-table slot. It reserves zero
for empty, uses a nonzero 30-bit state ordinal, a one-bit root selector and
a one-bit fingerprint. A matching fingerprint still triggers reconstruction
and comparison of the full canonical 83-bit root. Open addressing and the
first-inserted exact witness are preserved. The existing five-byte ID30 and
four-byte ID29 modes remain separately selectable.

For a state index `s` with `s+1 < 2^30`, selector `q < 2`, and fingerprint
bit `f < 2`, the stored word is

`f * 2^31 + (s+1) * 2 + q`.

It is nonzero and less than `2^32`. The low 31 bits decode `(s,q)` and the
high bit decodes `f`. `FourByteId30.lean` checks both implications with the
pinned Lean 4 toolchain. It also proves the 75% occupancy bound for the
`ceil(4*S/3)` request and at most two roots per state. The existing
`RootIndexCapacity.lean` checks the
N83 symmetric-state bound `K <= 5086` for the 30-bit state-ID mode. The
Rust unit test checks boundary encoding, first-witness retention, exact
lookup after slot/fingerprint collisions, and nonmatching keys. The bridge
from the Lean natural-number word to Rust `u32` is the checked bound; the
table's group-law behavior is checked by independent producer replay.

The next screen uses the previously sealed 128-point N83 batch. Freeze one
source and release binary, run K4096/root11 in fresh processes with an empty
base/index/log cache for both five-byte ID30 and four-byte ID30, and retain
the process order and exact flags. The four-byte table's expected allocation
at K4096 is 8 GiB versus 10 GiB for five-byte ID30; both have the same
state array, factor base and 128 public points. Set a 2,400-second cap for
each arm, retain timeout and failed replay, and independently reconstruct
the 679,936-point base, rank 4096, all 128 four-point relations and scalars,
target phase sums, memory and root-index layout. Compare per-point online
times and complete cold process separately. This ordinary-host screen is a
candidate diagnostic; a fresh-point Tier 1 same-point IC/rho panel follows
only if the exact four-byte mode yields a verified local gain.
