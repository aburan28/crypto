# The 30 B/s objective, priced

The objective carried through `README.md`, `RTX-PRO6000.md`, `COMPACT-STATE.md`,
`SHARED-SIGMA.md` and the `vector4-sigma` archive is *26 B complete scalar
iterations/s on one RTX PRO 6000*, raised here to 30. Six documents record it
as unachieved. None of them says what it would take, because the model that
could have said so — [THROUGHPUT-CEILING.md](THROUGHPUT-CEILING.md) — describes
the software multiplier that [native carryless](NATIVE-CARRYLESS.md) replaced,
and its verdict ("no configuration of this multiplier reaches 20 B/s") was
overtaken by the measured 14.637530 B/s.

This prices the objective against the path that actually runs.

## The boundary

The probe in [benchmarks/hardware-limits](benchmarks/hardware-limits/probe.cu)
measured what this GPU's integer pipe will issue: 62.1 lane-ops per SM-clock
for pure `LOP3`, 63.0 for `IADD3`, 69.3 for a `LOP3`/`IMAD` mix, and 78.3 for
`IADD3`/`FFMA` — the highest rate ever measured on the part. 188 SMs,
2.28–2.33 GHz sustained.

Put the optimisation history against that axis. Each row is a measured
throughput and a measured instruction-visit count from its own receipt:

| | B/s | instr/update | lanes/SM-clock |
|---|---:|---:|---:|
| software product, CUDA 13.0 preset | 7.110440 | 4,204.56 | **69.1** |
| + [native `clmad`](NATIVE-CARRYLESS.md), same layout, batch 32 | 8.703518 | 2,233.06 | **44.9** |
| + batch 16, [compact state](COMPACT-STATE.md), [shared masks](SHARED-SIGMA.md) | 14.637530 | 2,189.75 | **74.1** |

That tells the whole story of the last two rounds.

- The software path was **already issue-bound**: 69.1 against the measured
  69.3 for `LOP3` + `IMAD`, which is exactly the mix a software carry-less
  multiply produces.
- `clmad` cut instructions by **47%** and bought only **22.4%** of
  throughput, dropping the kernel to 44.9 lanes/SM-clock. It stopped being
  issue-bound; something else — state traffic and latency — became the
  binding constraint.
- The batching and layout work then bought **68%** at an essentially
  unchanged instruction count (×0.981). It did not remove work; it removed
  the stall, and put the kernel back on the pipe.

**The kernel is now at 74.1 lane-instructions per SM-clock: above every
pure-stream rate the probe measured and at 95% of the best mixed rate it ever
reached.** Both known mechanisms have been spent once each — `clmad` created
issue headroom, the layout work consumed it — and the kernel is back against
the integer pipe.

So from here throughput is a straight function of instructions per update, and
the objective converts exactly:

```
30 / 14.637530 = 2.05
```

**30 B/s on one RTX PRO 6000 requires cutting 2.05x of the instructions the
walk executes.** Halving every instruction in the kernel — all of them, with
no exception — lands at 29.3 B/s, still short. That statement needs no model
of where the instructions go; it follows from the measured rate alone.

The one caveat the table itself raises: a large instruction cut can move the
kernel off the pipe again, as `clmad` did, in which case the gain arrives only
after a second round of layout work. That is a reason to expect the next
factor to come in two steps, not a reason to expect more than the factor.

### A correction to the old model

THROUGHPUT-CEILING.md concluded that the kernel is integer-pipe bound "which
is why every instruction-count reduction in the history translated almost
one-for-one into speed and every scheduling, cache-policy, occupancy and
state-layout experiment did not."

The history since falsifies both halves. The largest instruction cut in the
project — `clmad`'s 47% — returned 22.4%, not one-for-one. And the largest
throughput gain in the project, ×1.682, came from exactly the
scheduling-and-state-layout experiments the model said do not work, at an
instruction count that did not move.

The model was not wrong about the mechanism; it was wrong to assume the kernel
is always on the pipe. It is on the pipe now, which is the only reason the
2.05x conversion at the top of this document is valid, and that has to be
rechecked after any change that removes a large block of instructions.

## Only ALU instructions count

The [vector4-sigma screen](benchmarks/vector4-sigma/README.md) is the
calibration. It replaced 56 scalar `LDS` with 14 `LDS.128`, cut instruction
visits from 2,189.75 to 2,147.75 — **1.9%** — and measured **+0.285%**
throughput, below its 0.5% qualification bar.

That is the expected result for a kernel saturating the integer pipe: the
instructions it removed issue on the LSU, which is not the resource that is
full. It also means the linear "each percent of instructions is a percent of
speed" rule in THROUGHPUT-CEILING.md holds only for ALU instructions, and
memory-instruction counts are close to free at the margin.

Everything below is therefore counted with ALU and memory apart.

## Where the update goes

`./path_cost.py` compiles the packed walk offline with clang and weights each
field routine by how often a scalar update runs it. The unit is clang PTX, not
the SASS the shipping nvcc build executes — `ptxas` fuses logic into `LOP3` and
the ratio is close to two — so the shares are the output, not the absolutes.

Batch 16, `clmad`, the audited RTX PRO 6000 preset:

| Component | per update | ALU/update | share |
|---|---:|---:|---:|
| inlined per-slot work (basis conversion, adds, squaring, weight, state, control) | 1 | 1,263 | 33.5% |
| polynomial product, paired (2 products) | 1.88 | 1,050 | 27.8% |
| **Frobenius network, walk (both coordinates)** | 1 | 695 | **18.4%** |
| normal-basis product, inverse chain | 0.50 | 336 | 8.9% |
| polynomial product, single | 1.06 | 307 | 8.1% |
| Frobenius networks, inverse chain | 0.31 | 122 | 3.2% |
| | | **3,774** | |

Memory instructions come to 129 per update against 3,774 ALU, which is why the
vector4 result looks the way it does.

Three things stand out.

**The products are 45% and they are the walk's floor.** 5.3125 field products
per update at batch 16 — 85 per batch, the five multiplications plus the
batched inversion the iteration function needs. That count is not negotiable
without changing the walk.

**The reduction costs more than the carry-less product.** Splitting one
polynomial product:

| | ALU | mem |
|---|---:|---:|
| `product131`, the 131x131 -> 262 carry-less product with `clmad` | 131 | 21 |
| `reducePolynomial131`, the direct reducer | **158** | 16 |

Native carryless made the multiply cheap enough that reducing the result is now
the larger half of a product. And that is not laziness in the reducer: **F(2^131)
has no irreducible trinomial** — checked exhaustively over all 130 middle
terms — so the modulus is dense and the generated reducer is a Barrett-style
fold rather than a two-XOR shift. The sparsest irreducible pentanomials are
`x^131 + x^8 + x^3 + x^2 + 1` and `x^131 + x^8 + x^5 + x^2 + 1`.

**A fifth of the update is the Frobenius network, and a further chunk is basis
conversion.** `sigmaWalkNetworkPairShared131` is 695 ALU for two coordinates,
and the forward loop runs four basis conversions per slot —
`fromPolynomial131` twice at 127 ALU and `toPolynomial131` twice at 159 —
another ~570 ALU, about 15% of the update. Neither is field arithmetic. Both
exist because the iteration function needs the **normal-basis** Hamming weight
while the products want a **polynomial** basis, and the walk crosses between
them twice per coordinate per step.

## What this says about the objective

30 B/s per GPU needs every instruction halved and then some. The products are
45% of the update and cannot be fewer, so halving the update means halving a
product too: 289 ALU down to about 145, which means the 158-ALU reduction has
to nearly vanish *and* the 131-ALU `clmad` product has to get cheaper. Against
a dense modulus forced by the field, that is not an engineering target.

**The per-GPU objective should be restated.** What the numbers support:

| Route | per-GPU rate needed | status |
|---|---:|---|
| 30 B/s on one GPU | 30.0 | needs > 2x on a kernel at 90–95% of measured peak issue |
| 30 B/s on two GPUs | **15.0** | **+2.5% from today** |
| 30 B/s on three GPUs | 10.0 | **already exceeded**; 3 x 14.64 = 43.9 B/s |

Two GPUs cross 30 B/s on a **2.5%** per-GPU gain. Three cross it today with no
code change at all, which [aws/README.md](aws/README.md) already budgets: "If
the goal is only the 15–20 B/s that one GPU could not reach, three GPUs do it."

## The levers, ranked by what is measured

Only ALU instructions pay, at roughly one percent of speed per percent removed.

1. **The reduction, 158 ALU of every product, 22% of the update.** It is the
   largest single addressable item and the one with a clear structural reason
   to be smaller. Changing to the `x^131 + x^8 + x^3 + x^2 + 1` basis would
   make the fold sparse; the cost is that the ONB conversion stops being
   structured, so this is only worth doing together with lever 2.
2. **The basis machinery, ~34% of the update between the Frobenius network
   (18.4%) and the four conversions per slot (~15%).** `fromPolynomial131`,
   `sigma^j` and `toPolynomial131` are all F2-linear, and the walk applies
   them in sequence. Their composition is one linear map, not three. Whether
   one general map beats the current structured three is exactly the
   measurement to make, and `path_cost.py` measures it without a GPU.
3. **The `clmad` product, 131 ALU.** `P131` holds 131 bits in five 32-bit
   words while `clmad` consumes 64-bit operands, so the product packs and
   unpacks around every call. Three 64-bit limbs would not.
4. **Not the memory path.** vector4 measured what that is worth: +0.285% for
   a 1.9% instruction cut.

## Reproducing

```sh
./path_cost.py                  # shipping preset
./path_cost.py --clmad 0        # the software-product path it replaced
./path_cost.py --target 30      # what an objective implies
```

It needs a CUDA-capable clang and the pip CUDA wheels that `ptx_stats2k.sh`
in `gpu/ecc2k` assembles a `CUDAPATH` from; no GPU, and no `ptxas` new enough
for `clmad`, because it counts PTX rather than assembling it.

## Checking the model against a known transition

The one place the composition model can be checked is the `clmad` step, which
has measured instruction visits on both sides. Running `path_cost.py` at the
batch that comparison used:

| | software | `clmad` | ratio |
|---|---:|---:|---:|
| this model, ALU PTX/update, batch 32 | 5,977 | 3,570 | 0.597 |
| measured SASS visits/update, batch 32 | 4,204.56 | 2,233.06 | 0.531 |

The model predicts the direction and most of the size of a 47% cut it was not
fitted to. It is 12% optimistic about how much `clmad` removed, which is the
expected sign: it counts a `clmad` as one instruction where the SASS count
sees the software multiply's `IMAD.WIDE` chain collapse further.

## What is not claimed

The unit here is clang PTX; the shipping build is nvcc SASS, and the two differ
by roughly a factor of two in absolute count. Shares are robust to that, ratios
between components are robust to it, and the boundary section does not use it
at all — it is derived from the measured rate and the measured instruction
visits in the receipts. The attribution of the inlined kernel body between the
two slot loops and the inversion chain is an estimate, and it moves the total
by about 8% either way; no conclusion above depends on it.

Lane-instructions per SM-clock are computed at the 2.30 GHz the probe
sustained. At the 2.43 GHz maximum the shipping row reads 70.2 instead of
74.1, which is still above every pure-stream rate measured.

None of the levers has been implemented or measured on a GPU here.
