# gpu/btcpuzzle — Pollard kangaroo for the interval ECDLP

A GPU solver for discrete logarithms **inside a known interval** on
secp256k1, and the Bitcoin puzzle registry that motivates it.

The Bitcoin puzzle addresses are a public challenge created in 2015: puzzle
#n holds a small amount at an address whose private key lies somewhere in
`[2^(n-1), 2^n)`. Knowing the interval is the whole game. Rho over the full
group costs 2^128 steps and is hopeless; Pollard kangaroo over an interval
of width W costs about 2·√W, so puzzle #n costs roughly 2^(n/2). Puzzles
into the 60s and 70s have been solved this way with open-source tools.

This reuses the verified secp256k1 field and point arithmetic from
`gpu/ecc`; what is new here is the kangaroo walk, the interval bookkeeping
and the solver.

**No kernel here has been run on a GPU.** The solver is verified end to end
on the CPU — it recovers planted keys through the full pipeline — and
`./bench selftest` checks the device path against that same host code.

## The two kinds of puzzle, and why only one is attackable

| What is known | Attack | Cost for puzzle #n |
|---|---|---|
| Public key (the output has been spent from) | Kangaroo, interval ECDLP | ~2^(n/2) |
| Address only (never spent) | Brute force the interval, hashing each candidate | ~2^n |

A Bitcoin address is HASH160 of the public key, so an unspent puzzle hides
its public key behind two hash functions. There is no interval algorithm
against a hash — you would have to walk the interval one key at a time,
squaring the exponent of the work. Every unsolved puzzle is far out of
reach that way.

**So this directory implements interval ECDLP only.** There is no address
scanner and no hash160 kernel, because for the puzzles it would not help,
and a bulk address scanner is a different tool with a different purpose.

## How the walk works

Shift the interval to be centred at zero: with `Q = k·G` and `k ∈ [a, a+W)`,
let `c = a + W/2`, `Q' = Q − c·G`.

```
tame kangaroo i    starts at  u_i·G,        distance u_i
wild kangaroo i    starts at  Q' + v_i·G,   distance v_i
```

Both herds take the same pseudorandom jumps: from a point P the index
`j = x_P mod NJ` picks a jump scalar `s_j`, and P advances to `P + s_j·G`
with the distance increased by `s_j`. When a tame and a wild kangaroo land
on the same point,

```
d_tame·G = Q' + d_wild·G    =>    k = c + d_tame − d_wild  (mod n)
```

Collisions are found through distinguished points, and kangaroos alternate
herd by index parity so every warp is half tame and half wild.

**Centring on the middle of the interval matters more than it looks.**
Both herds draw start offsets uniformly from `[0, W)`, so the tame herd is
centred on W/2. Measuring from `a` would leave the wild herd centred on
`k' + W/2`, so a key near the top of the interval would barely overlap the
tame herd at all. Measuring from the middle puts both herds in the same
place whatever the key is. Switching to the centre visibly improved the
measured constant at every interval size tried, most at the sizes where an
uncentred run had been worst; the before-and-after runs used different trial
counts, so the improvement is not quoted as a precise factor.

## The invariant that makes this testable

A tame kangaroo at distance `d` is at exactly `d·G`, and a wild one at
`Q' + d·G`. The tests check this **for every kangaroo after every step**,
which pins the jump table, the distance accumulation and the point
arithmetic simultaneously — a far stronger check than comparing against
recorded vectors, and it is also verified against the device's own state
in `./bench selftest`.

## Build and test

```bash
make test                       # CPU verification, no GPU required
make bench                      # CUDA solver (needs nvcc)
make bench ARCH=sm_100

./bench selftest                # kernels vs the host
./bench list                    # registry and cost model
./bench solve --bits 48 --self  # plant a 48-bit key and find it
./bench solve --bits 65 --pubkey 02...
./bench solve --puzzle 65       # target from puzzles.txt
```

## Measured cost

`make test` solves planted keys and reports the constant in
`steps = C·√W`. Reference points: 2.0 is the idealised two-herd bound, ~3.3
is Pollard's analysis of the classic lambda method, published parallel
implementations land near 2.1.

Kangaroo run lengths are close to exponentially distributed, so a mean over
a handful of solves is nearly meaningless — at 20 trials the standard error
is over 20% of the mean. `tune_kangaroo` reports the standard error
alongside the mean; on a 2^31 interval with 120 trials per configuration:

```
configuration                        mean  std err   median
apparent best (2^5, shift-2)         2.83     0.14     2.66
apparent worst (2^7, shift+0)        3.17     0.16     2.91
current default (2^6, shift+0)       2.79     0.13     2.65
```

So **C ≈ 2.8 ± 0.14**, between the idealised bound and the classic lambda
figure, and about a third above the best published implementations.

### Two things that were tried and did not work

**Tuning the jump parameters buys nothing.** A sweep over four jump-table
sizes and five mean-jump scales produced constants from 2.50 to 3.95, which
looks like a 36% win sitting there for the taking. It is not: at 20 trials
per cell that spread is entirely sampling noise, and re-running the apparent
best and the current default at 120 trials each puts them within half a
standard error of one another. The defaults are already on a flat optimum.
The `mean_shift` knob is kept so the measurement can be repeated at other
interval sizes, not because any setting of it is known to help.

This is the trap the table above exists to prevent: with twenty cells and a
20%-noisy statistic, the minimum of a grid search is a lucky draw, not a
tuning result.

**Restarting a kangaroo at a distinguished point is also a wash.**
Restarting unmerges kangaroos that collided within their own herd; not
restarting keeps the distance a kangaroo has built up. Across these sizes
the two are within noise, which is what a 64-entry jump table predicts —
same-herd merges are rare. Both are available; the solver defaults to
restarting.

## Occupancy: measured

`./ptx_stats_kangaroo.sh`, `k_kang_walk_lowmem<8>`, block size 128:

| Arch | `KG_MIN_BLOCKS` | Registers | Stack | Resident threads/SM |
|---|---|---|---|---|
| sm_90 | unset | 222 | 800 B | 256 |
| sm_100 | unset | 222 | 800 B | 256 |
| sm_120 | unset | 255 | 832 B | 256 |
| sm_90 / sm_100 | 4 | 128 | 992 B | 512 |
| sm_120 | 4 | 128 | 1328 B | 512 |

Constraining the launch doubles residency here, and it is worth more than
usual because `sm_120` otherwise pins to the 255-register ceiling — the same
consumer-Blackwell behaviour the prime-field rho kernels show, and the same
conclusion: tune `KG_MIN_BLOCKS` separately for consumer and datacenter
parts.

Seeding a kangaroo is what drives the stack frame, and it repaid the same
fix the rho kernels needed. Windowed scalar multiplication carries a
16-entry Jacobian table, 1.5 KB of frame that ptxas charges to every
resident thread for the entire kernel, to speed up an operation that runs
once per distinguished point. Swapping it for a table-free ladder:

| | Stack before | Stack after |
|---|---|---|
| `k_kang_walk_lowmem<8>` | 2320 B | 800 B |
| `k_kang_walk<8>` | 3264 B | 1792 B |

The jump table is staged in shared memory. A `kg_jump` is 25 words — eight
for the scalar, seventeen for the affine point — and 25 is odd, hence
invertible mod 32, so the entries a warp selects land in distinct banks.
That falls out of the layout rather than needing padding, but it is worth
knowing before anyone packs the struct.

## The registry

`puzzles.txt` holds one entry per line: the puzzle number, the compressed
public key (or `-`), and optionally a known private key. Every entry is
validated on load — the public key must decompress to a curve point, and a
supplied private key must lie in `[2^(n-1), 2^n)` **and** generate that
public key — so a bad paste is reported rather than silently searched for.

**The real puzzle public keys are not shipped here.** They are public
record, but this was written without a trustworthy offline source for them,
and a wrong key fails silently by searching for something that does not
exist. Paste them in yourself, or pass one with `--pubkey`. The entries
that *are* shipped are synthetic: keys generated locally, clearly labelled,
used by `make test` as regression targets. They are not puzzle keys.

## Files

| File | What it is |
|---|---|
| `kangaroo.cuh` | The walk: jump selection, distance accumulation, distinguished points, batched stepping |
| `kangaroo_host.hpp` | Jump table construction, SEC1 keys, the puzzle registry, collision → private key |
| `kernels_kangaroo.cuh` | Kernels and launch structure |
| `kangaroo.cu` | Solver and benchmark CLI |
| `test_kangaroo.cpp` | CPU verification, including end-to-end solves |
| `puzzles.txt` | The registry |
| `ptx_stats_kangaroo.sh` | Static instruction and occupancy analysis, no GPU needed |
| `tune_kangaroo.cpp` | Measures the constant with its standard error; slow, so not part of `make test` |

## Scope

This is an interval discrete-log solver. It does not scan address lists,
derive addresses, or touch a wallet or the Bitcoin network in any way — a
solved key is printed as a hex integer. The puzzle series is a deliberately
constructed public challenge, which is what makes it a reasonable
cryptanalytic target; the same solver applies to any interval ECDLP, such as
the biased-nonce and small-range key problems the rest of this repository
studies.
