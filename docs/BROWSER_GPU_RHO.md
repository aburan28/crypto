# Pollard rho in the browser, on the visitor's GPU

What the landing page runs when a visitor switches it on, what the number it
reports means, and what it deliberately does not do.

Files:

| Path | What it is |
|---|---|
| `docs/site/assets/rho-gpu.wgsl` | The walk as a WebGPU compute shader: 16-bit-limb Montgomery field, batched inverse, r-adding step, distinguishing test |
| `docs/site/assets/rho-gpu-host.js` | Host side in BigInt: jobs, curve arithmetic, jump table, walk seeds, trail replay, collision to logarithm |
| `docs/site/assets/rho-gpu-worker.js` | The device: buffers, dispatch pacing, the device-vs-host self-test, slot reseeding |
| `docs/site/assets/rho-gpu.js` | The page: the opt-in control and the readout |
| `scripts/rho_gpu_emulate.mjs` | Tests: the shader's arithmetic emulated in u32 JavaScript against BigInt, plus an end-to-end solve |

## 1. It is opt-in, and that is a design constraint

The control starts off, the choice is kept in `localStorage` for that browser
only, and a hidden tab pauses unless the visitor asked otherwise. Pausing
holds the live trails and the distinguished-point table, so nothing is lost
by it.

This is not only courtesy. A throughput figure taken from a visitor who did
not know their GPU was busy is not a figure worth quoting: nobody can say what
else the device was doing, and the measurement cannot be repeated by the
person it was taken from. `scripts/site/test_build.py` pins the default-off
state of both switches, because it is one HTML attribute away at all times.

## 2. The walk is the same walk

The step is the one in `gpu/ecc/rho.cuh` and `gpu/ecc/ecref.py`:

```
partition(P)  = limb0(internal(x)) & (R - 1)          R = 2^6 here
step          P <- P + M[partition(P)],  M[j] = c_j G + d_j Q
distinguished ((low32(internal(x)) >> 8) & dp_mask) == 0
```

`internal` is Montgomery form, as on the CUDA generic path. Two departures
from the CUDA engine, both deliberate:

- **No negation map.** It would buy a factor sqrt(2) in steps and cost the
  fruitless-cycle detection and canonical-element escape that
  `rho.cuh` needs to keep the walk a deterministic function of the point.
  This engine exists to measure a browser, so it measures the simpler walk.
- **16-bit limbs.** WGSL has no 64-bit integer type and no widening multiply,
  so a 32x32 product cannot be formed at all. With 16-bit limbs every CIOS
  term — accumulator word plus partial product plus carry — is at most
  `0xffff + 0xfffe0001 + 0xffff`, which is `0xffffffff` exactly. The bound is
  asserted on every term by the emulation test rather than argued in a
  comment.

### Occupancy is capped by the instance, not by the GPU

Every walk is abandoned at its first distinguished point and only the host can
reseed it, so the trails alive at any moment carry a combined tail of
`slots * 2^dp_bits` steps that no collision is ever found in. Filling a large
GPU with a small instance spends most of its work there: 8192 walks on toy32
overshot the expected `2^15.7` steps by 8.7x before the occupancy was capped
at `2^(rho_log2 - dp_bits)`. With the cap, toy32 runs 512 walks and finishes
around 1.9x its expectation, which is the ordinary distinguished-point
overhead rather than an artefact of the occupancy. The per-dispatch step count
is capped for the same reason, at half the distinguishing interval.

One inversion per `BATCH = 8` walks: the denominators of eight concurrent
steps are inverted together by Montgomery's trick, so Fermat inversion
(`a^(p-2)`, the only inversion WGSL makes easy) is amortised eight ways. Each
thread therefore owns eight walks, and 8192 walks is the default occupancy.

## 3. A distinguished point is four words

The device reports `(slot, restart, steps, x)` and keeps coordinates
otherwise. A trail's start is `aG + bQ` with `(a, b)` drawn from the job seed
by `(slot, restart)`, so the host can rebuild any trail — and its
coefficients — from those two indices, and does so only for the two trails
that actually collide, at a cost of about `2^dp_bits` BigInt point additions
each. That is the same argument the CUDA host makes, and it is why `dp_bits`
is set per instance in `rho-gpu-host.js`: large enough that the
distinguished-point table stays small, small enough that a replay is seconds
rather than hours.

A collision on `x` is a collision up to sign, so both signs are tried and the
logarithm is reported only if `kG = Q` verifies.

## 4. What the measurement is, and is not

The instances are toy curves from `gpu/ecc/ecref.py` (`make_toy_curve`, seed
7) at 32, 40, 48, 56 and 64 bits, with prime group order and a planted `Q`
whose logarithm is not shipped. Expected walk length is `sqrt(pi n / 4)`,
which the page shows as the exponent the run is working against.

So the figure this produces is **throughput**: steps per second, and the
expected number of steps it has to be divided into. It is not a ratio to the
generic-group floor and not an `S = ops / sqrt(n)` entry in the cost ledger —
the curves here are far below any size where that comparison says anything,
and the reason to run them is that a browser can finish one. A solver whose
inner loop is fast and whose collisions do not solve is a common way to
publish a throughput number that means nothing; running these to completion
is what rules it out.

Nothing here approaches a deployed curve. A 64-bit instance is about `2^31.4`
steps; secp256k1 is about `2^128`.

## 5. What is tested, and what is not

`node scripts/rho_gpu_emulate.mjs` re-implements the shader's arithmetic in
u32 JavaScript, function by function, and checks it against BigInt on every
job curve: Montgomery multiplication and its carry bounds, the conditional
add and subtract, Fermat inversion, the batched inverse, the partition and
distinguishing predicates, and 64 consecutive steps of the real walk against
a BigInt walk on the real jump table. It then solves toy32 end to end through
the emulated device and the engine's own replay path — about 10^5 steps, one
collision, `k` verified.

**It does not compile WGSL.** Nothing in this repository can claim that
`rho-gpu.wgsl` matches that emulation, in the same way that `gpu/ecc` cannot
claim its kernels have run on a GPU. The gap is closed at runtime instead:
before the worker posts a single rate, it replays its first 32 device steps
for 16 walks in BigInt on the host and refuses to run if the device disagrees.
The page shows the outcome of that self-test next to the device name, and the
build test pins the ordering so a future edit cannot report a rate first and
check the device afterwards.

## 6. What has actually been run

The engine was run end to end in headless Chromium against Chromium's
**software** WebGPU adapter (`--use-webgpu-adapter=swiftshader`), which is the
only adapter available in the container this was developed in:

| Instance | Walks | Steps to the solve | Wall clock | Logarithm |
|---|---|---|---|---|
| toy32 | 512 | 108,892 (1.9x the expectation) | 1 s | recovered, `kG = Q` verified |
| toy40 | 1024 | 2,471,248 (3.3x) | 10 s | recovered, `kG = Q` verified |

Both recovered the `k` that `ecref.py` planted. That is what the table
establishes: the shader compiles, the device agrees with the host, and a
collision between two device trails yields the logarithm.

What it does not establish is a GPU throughput figure. SwiftShader is a CPU
rasteriser; roughly 250k steps/s there says nothing about real hardware, and
no rate from it is published on the site. The rate a visitor sees is measured
on their own device, and is the only rate this engine claims.
