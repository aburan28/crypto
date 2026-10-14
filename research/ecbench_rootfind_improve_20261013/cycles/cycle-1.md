# Cycle 1 — |F| select + rootfind-window

Session: `sessions/select` (`ECBS1h5cb3bc5c9cd5`), 42/42 verified, audit OK,
replay cert `8421ef064ead9419c7ab8a520ac8f968fdd0f9fe37c75e1ef4474f249460401a`.

## Mean cold S (n=13, 2 targets × 2 rounds)

| arm | mean S | setup GAE | online GAE | \|F\| | solve wall (ms, L0) |
|---|---:|---:|---:|---:|---:|
| rho-strong | 14.0 | 0 | — | — | 0.10 |
| ic-f4-m3-d5 | 15.6 | 0 | 695 | 33 | **191** (lower-bound S) |
| ic-rootfind-m2-d4 | **30.2** | 190 | 1161 | 19 | 0.56 |
| ic-rootfind-m2-d5 | 31.6 | 561 | 849 | 33 | 0.78 |
| ic-rootfind-m2-d3 | 32.6 | 66 | 1394 | 11 | 0.72 |
| ic-rootfind-window-m3-d5 | 32.9 | 561 | 908 | 33 | 0.65 |
| ic-rootfind-m3-d5 | 39.2 | 561 | 1191 | 33 | 0.84 |

## Decisions

1. **Selected cold candidate:** `m=2`, `dimension=4` (lowest cold S among
   rootfind arms; setup 190 vs 561 at d=5).
2. **Keep window engine:** `rootfind-window` at d=5 cuts mean S 39.2 → 32.9
   vs full m=3 (same table build). Not stream-identical to F4.
3. **d=3 is a dead end for cold S** on this smoke: tiny setup, but online
   trials erase the gain.
4. **Amortization (secondary):** sharing setup across N=8 helps d=5 most
   (31.6 → ~20.6 equivalent); still above rho. Labelled secondary only.
5. **Wall:** every rootfind arm ≪ F4 m=3 descriptively (~0.6 ms vs ~191 ms);
   L0, not admitted. Do not sell as S speedup.

## Next (cycle 2)

Measure `rootfind-window` at d=4 (smaller table + window). Run the claim
panel with d4/m2 + window arms vs F4 m=3 and #1765 d5 baselines.
