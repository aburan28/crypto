# N83 compact exact-root table at the K1102 capacity edge

The opt-in `KIC_INDEXED5_COMPACT27=1` requests exactly 2^27 five-byte
slots for the symmetric one-root N83 index. At K1152 this is 640 MiB,
versus 1,280 MiB for the existing dynamic load75 request. Every occupied
slot still stores an exact state/root witness with the same fingerprint and
key reconstruction. The producer refuses this mode unless the complete
two-root bound is at most 85% of the slots, the layout is the one-root
five-byte index, and the Bloom option is off. The default layout and its
capacity rule remain unchanged.

The initial experiment uses K1152, rank seed 29, the already exposed
one-point N83 public input at SHA-256
`2c068d2f659588edcdbaa685fc73bb146e2bfee3be123d377d537658f8021c6d`,
and exact IC source SHA-256
`e1b143ec8ba80425a1102b65c3e7c0b67a4b4e9a7e0bc116bcaf13e76c5651f9`.
The same release binary runs four fresh processes in mirrored order:
`off_a`, `compact_a`, `compact_b`, `off_b`. All other source, target,
factor-base, root-index and rank flags are identical. Each child has a
900-second cap. Preserve stdout, stderr, full enumerated base, recovered
target, exact command and environment, outer wall, macOS peak RSS, and
producer setup/rank/target/process intervals.

`verify.py` independently reconstructs every base, sums the four factor-base
points for the public target, replays the recovered scalar, checks the timing
phase boundary, and requires the same base hash and exact root-entry count
across all arms. It reports the two A/A spreads and full cold-process means.
This local host supplies an algorithm and correctness diagnostic. Advance
the compact mode to a Tier 1 paired cloud comparison only if all four arms
replay, both A/A spreads are at most 10%, and the compact mean cold advantage
exceeds the largest absolute A/A difference. Preserve a failed or noisy
panel with its exact status and cap.
