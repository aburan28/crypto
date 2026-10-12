# Protocol E-1: a packed root index for the compact-orbit producer

Frozen 2026-10-07, before any instrument is built.  **Engineering**, by
§3 of `AGENTS.md`: it moves bytes and probes, not an exponent.  `S` is
reported before and after; no speedup against rho is claimed.  Status:
**PENDING**.

## Derivation (stated before measuring)

The compact-orbit producer keeps one `S₃` root index of `K²·n` states and
answers each probe with one lookup.  Probes per relation scale as
`r / (n²K²)`, so at fixed memory the index's bytes per state set `K`, and
`K` sets the probe count in both the rank and the online stages: the
repository's own law is five bits of time per ten bits of memory.

The current layout, read from the unmerged wide-path producer on branch
`cursor/ic-boundary-experiments-d111` (`examples/koblitz_orbit_dlp_fast.rs`,
commit `7b1787074`, labelled unmerged):

| structure | per entry | why |
|:--|--:|:--|
| `RootTable128.slots: Vec<(u128, u64)>` | 32 bytes, padded | open addressing, slot count the power of two at or above `2 × entries`, so at least 64 bytes per stored entry |
| `State128` | about 40 bytes | three `u16` indices plus two `u128` normal-basis roots |
| transient `shifted` vectors and labels | unmeasured | |

The one measured total is the `n = 83`, `a = 1` run: 29,880,000 regular
states inside a 5,771,837,440-byte peak, about 193 bytes per state.  The
root table and states alone account for at least 104; the rest is
unattributed and is measured first.

**Proposal.**  Replace the open-addressed table by a sorted array of packed
entries: the canonical root in `n` bits and the state id in
`⌈log₂ states⌉` bits, in one `u128` whenever `n + ⌈log₂ states⌉ ≤ 128`
(true to `n = 97` with `2^{31}` states; `n = 127` needs 24 bytes), with a
radix directory on the top 20 bits of the root (`2^{20}` offsets, 4 MiB)
so a lookup is one directory read plus a short bounded scan.  Drop the
stored roots from `State128`: on the rare hit they are recomputed from the
three indices.  Expected: 16 bytes per entry plus 6 bytes per state, under
32 bytes per state including the directory.

## Instrument (Rust, to build in the follow-on PR)

1. **Attribution first.**  Instrument the producer to report resident
   bytes per index structure at the end of construction, and peak RSS, at
   `K = 600` on `icv1-f2m83-t6151469093347-cdcc5432` and at `K = 440` on
   `icv1-f2m53-tm56619371-dac20a85`.
2. Implement the packed index behind a flag (`KIC_INDEX=packed|hash`),
   same scan order, same canonicalisation, same first-hit rule.
3. Run both index kinds at equal `K` on both curves, three runs each,
   under `tools/isolated_bench.py` where the host allows it; otherwise
   report instruction counts and mark wall time contended.
4. Grow `K` under the packed index until resident bytes reach the hash
   index's envelope at `K = 600`, and run the rank stage and the frozen
   `R1` online target at that `K`.

## Predictions (pass/fail)

- **M1 (bytes).**  Packed index resident bytes per state at most 32 at
  `n = 83`, against at least 104 for the hash layout.
- **M2 (rate).**  Probe rate within `[0.75, 1.25]` of the hash index at
  equal `K`, in instructions per probe.
- **M3 (regression).**  At `K = 600` the packed index reproduces the
  frozen `R1` relation byte for byte, including the 8,845,441-probe count.
- **M4 (K at equal memory).**  Under the hash index's envelope the packed
  index runs `K ≥ 1,200` at `n = 83`.
- **M5 (probes).**  At that `K`, rank-stage mean probes per relation fall
  by at least 3.5× against `K = 600`, and the online probes to the same
  frozen target fall by at least 3× or the first-hit rank is explained by
  the scan order.

## Decision rule (registered)

M1, M2 and M3 passing makes the packed index the default; M4 and M5 size
the gain.  The row lands on the frontier with its memory axis and its
instruction count beside the strong rho, classed **engineering**.  Any
failure of M3 is a bug, not a result.

## Stop condition and inadmissible moves

Bounded: two curves, two index kinds, one `K` sweep.

Inadmissible: changing the scan order or the first-hit rule while packing;
reporting the probe reduction without the memory it bought; comparing the
new `K` against the old rho envelope rather than the same one; reading any
of this as a change to the ratio-to-rho exponent, which the strong-rho
sweep measured at `+0.078` per doubling of `r` and which this does not
touch.
