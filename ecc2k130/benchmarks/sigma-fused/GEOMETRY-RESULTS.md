# Fused sigma geometry: retain batch 16

The fused schedule does not make a larger batch competitive on the RTX PRO
6000. All five arms passed 300/300 host replay, zero drops and identical sorted
1,710,327-record v1 corpora. Every B32/B64 timing sample was slower than every
B16 sample; no geometry qualifies for confirmation.

Each screen sample completes the same 201,863,462,912 scalar updates from
6,160,384 live walks. Five warmups are excluded. Three rounds used the frozen
forward/reverse arm orders. This is a screen, with a separate GPU allocation
from the confirmed 15.436677 B/s B16 comparison.

| arm | median B/s | median / B16 | range B/s | replay/corpus | decision |
|---|---:|---:|---|---|---|
| B16/T256/min2 | **15.275127** | 1.000000 | 15.249411–15.384403 | pass | retain |
| B32/T256/min2 | 13.695901 | 0.896615 | 13.693160–13.697531 | pass | reject |
| B32/T512/min1 | 13.458151 | 0.881050 | 13.453548–13.459296 | pass | reject |
| B64/T256/min2 | 10.663847 | 0.698118 | 10.660661–10.665862 | pass | reject |
| B64/T512/min1 | 10.640139 | 0.696566 | 10.628282–10.650520 | pass | reject |

The native post-run checker reopens all twenty warmup/screen logs, validates
their arm, launch bounds, full backend markers, final counted work, rates and
zero drops, then scans the retained sorted corpus for all 1,710,327 records,
run-id 31 and canonical 131-bit keys. All twenty log hashes and eight measured
source hashes match. Recompiling the native summarizer reproduces `result.json`
byte for byte. These checks establish the bounded rejection; they do not add
a gain to the earlier confirmation.

The first attempt remains preserved as a producer failure with no timing: its
replay and corpus checks passed, but it required an automatic-worker occupancy
line that explicit-worker runs do not print. The additive repair introduces a
compile-time launch-bounds marker and retains the same populations and protocol.

The final artifact and Modal retrieval details are in
`results/geometry/artifact.json`; its archive SHA-256 is
`785fa37ac4ba2df15f0d36dbd6638ab6712bd5b09d7813c6c126d5c50b21f60b`.
No arm changes the sigma walk. The active 26 B/s objective remains unmet.
