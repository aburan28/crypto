# Iteration function and the 28 B/s question

Question: can a different iteration function take one RTX PRO 6000 from the
audited 14.64 B updates/s to 28 B/s?

Answer: **no.** 28 B/s is below the floor of *any* walk that performs one
affine point addition per iteration, which is every walk that keeps the rho
collision structure. The floor is set by the carry-less unit, not by the
Frobenius work the iteration function controls. The redesigned function of §4
removes the Frobenius and basis-conversion work; built and measured (§6) it
cuts the static ALU cost per update by 21% and runs **16.4% more updates per
SM-clock** than the shipping walk on an RTX PRO 4500 (8.7% more at that card's
165 W power limit), with every device report re-walked correctly. Carried to
the 6000 that is 15.9–17.0 B/s, an extrapolation that straddles the 17.0 B/s
target of §4.4 and stays undecided until a 6000 is available. The rest of the
distance to 28 is field products, and products do not move with the walk.

This note states the boundary before anything is optimised (`AGENTS.md` §1),
prices the current kernel and the floor in one unit (§2), records the accounting
corrections that came with it (§3), specifies the replacement walk with its
falsification target (§4), and reports the measurement against that target
(§6).

## 1. Boundary

Unit: **lane-instructions per scalar update**, split by execution pipe, on the
188-SM RTX PRO 6000. Both pipe rates are measured, not quoted:

| pipe | rate | source |
|---|---:|---|
| integer/logic (`LOP3`, `IADD3`, `SHF`, `PRMT`, …) | 62.1 lanes / SM-clock | `benchmarks/hardware-limits/probe.cu` on the 6000, [THROUGHPUT-30B](THROUGHPUT-30B.md) |
| carry-less (`CLMAD`) | 1 / 37.7 of the logic rate = **1.65 lanes / SM-clock** | `benchmarks/clmad-price/probe.cu`, measured on an RTX PRO 4500 (same SM) |
| sustained clock | 2.30 GHz (probe), 2.43 GHz maximum | THROUGHPUT-30B |

Per-update budgets follow from `SM-clocks/update = 188 × f / rate`:

| target | SM-clocks / update | ALU budget (100%) | clmad budget (100%) |
|---:|---:|---:|---:|
| 14.64 B/s (measured) | 29.5 – 31.2 | 1,834 – 1,938 | 48.7 – 51.4 |
| 20 B/s | 21.6 – 22.8 | 1,343 – 1,418 | 35.6 – 37.6 |
| **28 B/s** | **15.4 – 16.3** | **959 – 1,013** | **25.4 – 26.9** |

(ranges are 2.30 → 2.43 GHz; budgets are at 100% pipe utilisation, which no
kernel reaches — 90% is the practical line.)

The clmad-rate row is the one number here not measured on the 6000 itself; it
is the first thing to re-measure when a 6000 is available, and every clmad
ceiling below scales with it.

### The floor

An iteration of any rho walk on this curve that keeps the class structure is
one affine addition `R + Q` with a per-lane `Q`, and affine addition with
batched inversion costs a fixed number of field multiplications:
`λ = e/d` (1), `λ(x + x')` (1), the Montgomery chain (3, asymptotically), plus
one squaring. At batch 16 the kernel issues **77 polynomial products per 16
slots = 4.81 per update**, and a 131×131 carry-less product by 2-way Karatsuba
is **6 `CLMAD`** (three 64×64 products, `lo` and `hi` each). Nothing in the
iteration function changes either number: `Q` can be `σʲ(R)`, a table point, or
anything else, and the addition still costs the same.

| component | clmad / update | ALU / update (static SASS) |
|---|---:|---:|
| 4.81 products (Karatsuba + dense reduction) | 28.9 | ≈ 830 |
| squaring in the polynomial basis | 5.0 | ≈ 80 |
| batched inversion, amortised over 16 | 4.6 | ≈ 270 |
| **arithmetic floor** | **38.5** (28.9 with squarings moved to ALU) | **≈ 1,090 – 1,180** |

Set against the budgets:

- **clmad**: 28.9 per update is the absolute floor (squarings and the inverse
  chain's squarings pushed onto the ALU). At 1.65 lanes/SM-clock that is
  17.5 SM-clocks, i.e. a ceiling of **24.7 – 26.1 B/s at 100% utilisation of
  the carry-less unit with every other cost free.** 28 B/s needs ≤ 25.4 – 26.9
  clmad per update *including* squarings and inversion. There is no room.
- **ALU**: the arithmetic alone is ≈ 1,090 static lane-instructions, already
  above the 959 – 1,013 budget at 28 B/s, before the hash, the canonicalisation
  that makes the walk class-invariant, the table or Frobenius work, state
  traffic and control are paid.

Two independent pipes, two independent floors, both above 28 B/s. This is the
same shape as the residual-walk thread's §11.7: a target quoted from one phase
(the Frobenius/conversion share the walk controls) that the unpriced phase
(products) forbids.

What would move the floor, and why each is out of reach here:

| lever | effect | status |
|---|---|---|
| fewer than ~4.8 products per step | lowers both floors | no affine formula does it; projective forms need the affine `x` for the hash anyway (an inversion), λ-coordinates cost 1I+3M+2S, `2P+Q` single-inversion tricks cost ≥ 8M |
| a product in fewer than 6 clmad | lowers the clmad floor | 32-bit limbs need 9 `lo` products; 3-limb splits need `hi` anyway; no 128-bit `CLMAD` form exists |
| products off the carry-less unit | trades pipes | software `clmul` costs ~9 SM-clocks per product on the FMA pipe vs 3.9 on clmad; bit-sliced ALU products cost ~300 ALU each and the ALU is the tighter pipe; tensor cores need a shared matrix operand, and every walk's operands differ |
| a walk that is not one addition per step | changes the count | `x`-only doubling is 1I+1S but a single multiplier `[2]` is a permutation, not a random function; mixed doubling/addition diverges per lane, and SIMT pays both |

## 2. The single table

One unit, every variant a row, the floors as the last rows. Columns are static
SASS lane-instructions per update from `kernel_cost.py` (§3), with quarter-rate
`POPC`/`FLO` priced at 3.97 slots and `IMAD` at 2.01 (§3, measured). Measured
rates are RTX PRO 4500 medians from
[`benchmarks/table-walk/comparison.json`](benchmarks/table-walk/comparison.json),
reported both at the card's power limit and per SM-clock; the 6000 column is
an **extrapolation** by the per-SM-clock ratio (§6.3) and is marked as such.
Superseded figures stay as "before" marks.

| variant | ALU slots / update (static) | clmad / update (static) | 4500 B/s (median of 3) | updates / SM-clock | 6000 B/s | ALU / floor | correct | class |
|---|---:|---:|---:|---:|---:|---:|---|---|
| shipping walk `R ← σʲ(R) + R`, batch 16, audited preset | 2,271 → **2,324** (quarter-rate priced) | 45.25 | **5.107** @ 1.95 GHz | 0.0319 | **14.64** measured | 2.1 | 300/300 | baseline |
| table walk, bit-plane selection (first build, superseded) | 2,147 | 45.25 | 5.493 | — (clock not sampled) | — | 2.0 | 300/300 | engineering |
| **table walk of §4, byte/nibble LUT selection** `R ← R + ε·σᵏ(T_h)` | **1,831** | 45.25 | **5.552** @ 1.82 GHz (+8.7%) | **0.0372** (+16.4%) | *15.9 – 17.0 (extrapolated)* | **1.7** | 300/300, 64/64 planted logs | **engineering** |
| table walk, LUT selection + squarings on ALU | 1,938 | 40.25 | 5.522 (+8.1%) | 0.0368 | — | 1.8 | 300/300 | engineering (neutral: ALU is the binding pipe) |
| arithmetic floor, one affine addition per step | ≈ 1,090 | 28.9 | | | 22.2 – 23.4 (90% of either pipe) | 1.0 | | floor |
| **28 B/s budget** | 959 – 1,013 | 25.4 – 26.9 | | | | 0.88 | | below floor |

Reading the table: the ratio column is ALU per update over the arithmetic
floor. The shipping walk is at 2.1× the floor; the walk change measured takes
it to 1.7×; the predicted 1.3× (previous revision of this table:
"≈ 1,400 – 1,500") was not reached because the walk's own selection logic
costs ~175 slots and the arithmetic's share was under-priced before quarter-rate
instructions were measured (§3). 28 B/s is at 0.88× — a ratio that no
iteration function can produce because the floor is not made of
iteration-function work.

Per `AGENTS.md` §3, every measured walk row is **engineering**: the ALU/floor
column fell from 2.1 to 1.7, but the floor here is the cost of one affine
addition, i.e. of the *same* generic algorithm, so the column measures
overhead removed and not a generic-group bound crossed. The expected number
of iterations is unchanged up to the r-adding constant (§6.2), and the ratio
is bounded below by 1.

## 3. Accounting correction

`sass_cost.py` prices routines in isolation and, on the shipping preset,
reported the kernel at ~1,570 ALU and 29 clmad per update. Both were low:

- the audited build compiles the product and Frobenius routines as
  `__noinline__` (`ECC_BIG`), so the kernel body contains `CALL`s and the
  routine bodies are counted zero times by a body-only count;
- the paired product (`mulPolynomialPair131`) is two products and was weighted
  as one.

`kernel_cost.py` (added with this note) compiles `walk()` for `sm_120` with the
preset, dumps SASS, finds the step and slot loops from backward branches,
attributes each callee's body to its call sites, and divides the inversion and
step overhead by the batch. On the audited preset it reports **2,271 ALU-slot
and 45 clmad static lane-instructions per update**, 2,546 instructions in all.

The static figure exceeds the ALU pipe's capacity at the measured 14.64 B/s
(1,834 – 1,938 lane-instructions per update at 100%), so the dynamic count is
at most that: predicated-off arms of the slot-0/pair branches account for the
difference. Either way the kernel is at the ALU pipe's limit, which is what
Nsight's 84% ALU utilisation on the 4500 said and what THROUGHPUT-30B concluded
from PTX shares. The clmad pipe is at 40–45 per update against a 48.7 – 51.4
budget: **75 – 90%**. Both pipes are near saturation at 14.64 B/s, which is
why THROUGHPUT-30B's "clmad is one-for-one" expectation did not hold and why
any change that trades ALU for clmad (the `spread32p` route) buys nothing.

Class: **accounting.** Numbers changed, the kernel did not, and no gain is
claimed.

### 3.1 Quarter-rate instructions

The first table-walk build (bit-plane selection, §6.1) was predicted at
≈ 1,900 static ALU and measured only 7.6% faster than the shipping walk; the
gap was in the pricing of `POPC` and `FLO`, which `kernel_cost.py` had counted
as one slot each. Measured on the 4500 with the extended
`benchmarks/clmad-price/probe.cu` (log in
[`benchmarks/table-walk/raw/rtx4500-instruction-rates.log`](benchmarks/table-walk/raw/rtx4500-instruction-rates.log)),
lanes per SM-clock at 2.407 GHz on 82 SMs:

| instruction | lanes / SM-clock | slots (LOP3 = 1) |
|---|---:|---:|
| `LOP3`, `SHF`, `PRMT` | 64.0 | 1.00 |
| `ISETP` + `SEL` pair | 95.0 (pair counted as two) | 0.67 |
| `IMNMX` | 35.4 | 1.81 |
| `IMAD.WIDE` | 31.6 | 2.03 |
| **`POPC`, `FLO`** | **16.0** | **3.99** |
| `LDS.U8`, random shared address | 9.2 (MIO pipe, not ALU) | — |
| `CLMAD.lo` | 1.687 | 37.9 |
| 64×64 carry-less product (`lo`+`hi`) | 2.00 products | 63.9 |

`kernel_cost.py` now prices `POPC`/`FLO`/`BREV` at 3.97 slots. The shipping
walk carries 17.8 of them per update (weight, distinguished-point test, hash),
so its static figure moves 2,271 → **2,324**; the bit-plane selection carried
72 (its masked popcounts), which is why its 1,934 nominal was 2,147 real and
why the LUT selection of §6.1 replaced it. Class: **accounting**, and the
correction is what made the second design round land.

## 4. The replacement walk

What the iteration function *does* control is the third of the update that is
not arithmetic: the Frobenius network on both coordinates (`σʲ(x)`, `σʲ(y)`,
~460 ALU), two `fromPolynomial131` (~120) and two `toPolynomial131` (~160), in
all ≈ 740 static ALU per update, existing only because the shipping walk needs
the **normal-basis** weight of `x` for the branch and then applies `σʲ` to a
point whose products want the **polynomial** basis.

### 4.1 Definition

Walk on the classes of `⟨σ, −1⟩` (262 points each) as today, with a table of
`H` points `T_h = a_h P + b_h Q` of known coefficients:

```
h(R)  = (HW(x_n) / 2) mod H                      normal-basis weight, class-invariant
k(R)  = (Σ_e L(e) · x_e) · HW(x_n)⁻¹  mod 131    Frobenius phase
ε(R)  = bit p(R) of y_n, p(R) = argmax_{e ∈ supp(x_n)} (L(e) − k(R)) mod 131
R'    = R + (−1)^ε · σ^{k(R)} (T_{h(R)})
```

`L(e)` is the discrete logarithm base 2 of the folded coordinate index `e`
in `(Z/263)*/±1`, a fixed 131-entry table; `x_n`, `y_n` are the coordinates
in the permuted normal basis. The three identities that make this a class
function:

- `k(σR) = k(R) + 1 (mod 131)`: squaring sends coordinate `e` to `fold(2e)`,
  so `L` shifts by one at every set bit and the sum shifts by `HW(x)`;
  dividing by `HW(x)` (nonzero mod 131 since `0 < HW(x) < 131` on the curve)
  makes the shift exactly one.
- `ε(σx, σy) = ε(x, y)`: `p` is the support element "last before the phase" in
  the cyclic `L`-order, which is carried along by σ.
- `ε(x, y + x) = 1 − ε(x, y)`: `p` is a set bit of `x`, so negation flips it.

Hence `f(σR) = σR + (−1)^ε σ^{k+1} T_h = σ f(R)` and
`f(−R) = −R − (−1)^ε σ^k T_h = −f(R)`: the walk descends to classes with no
fruitless-cycle handling, exactly as the shipping walk does. This is an
r-adding walk (Teske) on the class set; `H = 16` gives the mixing the 8-branch
Frobenius walk has today, and the expected-iterations figure `2^60.9` in
`README.md` moves only by the r-adding constant, which the small-curve harness
measures directly.

All three identities, plus `Tr(ab) = parity(a ∧ b)` in the permuted ONB (the
trace form is self-dual, which makes trace tests of products free), were
verified on the code generator's `Onb` model at `m = 23, 41, 83, 131` over 300
random `(x, y)` each with zero failures.

**Why negation is exact and not probabilistic.** Two trails that meet as `R`
and `σᵏ(−R)` must stay in one class for the ~2^25 steps to the next
distinguished point; a canonicalisation that fails on a fraction `p` of classes
splits them after `1/p` steps, so `p` must be far below 2^−25. Weight
comparisons and other cheap invariants tie at rates of 10% and cannot be used;
`ε` above is complete. The algebraic alternative `Tr(y/x)` is also complete but
needs `1/x`, which the batched inversion of `x + x_T` does not provide.

### 4.2 Solver and corpus

The endpoint of a trail is `(α₀ + Σ_t ε_t s^{k_t} a_{h_t}) P + (β₀ + …) Q`;
the exponents are recovered by replaying the trail from its seed, as now,
with the table coefficients in place of the product of `(1 + s^j)`. The
distinguished-point predicate (`HW(x_n) ≤ 34`) is unchanged and stays
class-invariant. **The walk is a different function, so distinguished points
from the shipping walk and from this one cannot collide.** `docs/ecc2k130-status/history.json`
records zero collected points at the time of writing, so the fork costs
nothing today; once collection starts under either function it is fixed for
the campaign.

### 4.3 Cost

Per update, static SASS estimates against the shipping preset:

| removed | ALU | added | ALU |
|---|---:|---|---:|
| Frobenius network, two coordinates (695 PTX-ALU in THROUGHPUT-30B; SASS fuses at ≈ 1.66:1) | −420 | `k(R)`: eight `L`-bit-plane masked popcounts over 5 words, fold, mod 131, inverse lookup | +100 – 120 |
| `fromPolynomial131(y)` (127 PTX-ALU) | −75 | `p(R)`: threshold mask `{e : L(e) < k}` (131-entry table) and an 8-plane greedy argmax over the support, one `LOP3`/`POPC` pair per plane per word | +100 – 120 |
| `toPolynomial131(d)`, `toPolynomial131(e)` (159 PTX-ALU each; `d = x_p + x_T,p`, `e = y_p + y_T,p` are formed in the polynomial basis directly) | −190 | `fromPolynomial131(y)` for the `ε` bit, `ε`-select, table adds | +95 |
| shared σ mask loads | −(mem) | table lookup: 3 × `LDS.128`/`LDS.32` per slot from 131 × H × 40 B (84 KB at `H = 16`, 42 KB at `H = 8`) | +(mem, ~10 ALU) |
| | **−685** | | **+305 – 345** |

Net **≈ −340 to −380 static ALU**: 2,271 → ≈ 1,900 static, ≈ 1,450 – 1,550
dynamic by the same static→dynamic ratio as today. At the ALU pipe's rate that
is **16 – 18 B/s at 90% utilisation**, with clmad at 40 (17.1 B/s ceiling at
90%) unless the polynomial squaring is moved to the ALU (+50 ALU, clmad 35,
ceiling 19.6). The honest range is **16 – 19 B/s**, and the clmad row is the
binding one at the top of it.

*Outcome (§6):* the estimate above was written before quarter-rate pricing
(§3.1). Measured, the shipping walk is 2,324 slots and the built table walk
1,831, a net **−493 slots** — more than the −340 to −380 predicted, because the
LUT selection replaced the popcount-based `k`/`p` of the "added" column with
shared-memory lookups (+≈ 175 slots, and 126 memory instructions per update
against the shipping walk's 107 including its shared σ-mask loads, instead of
+305 – 345 slots). The clmad count did not change (45.25 static; the squaring stays
on the carry-less unit), and moving it to the ALU measured neutral because the
ALU, not clmad, is the binding pipe at the resulting rate.

### 4.4 Falsification target

The redesign is worth landing iff, on the RTX PRO 6000 with the audited
geometry and the `dpWeight = 34` predicate:

- benchmark median **≥ 17.0 B/s** (a 16% gain, the point where the fork's
  cost is clearly repaid) — below that, classify as **engineering that did not
  pay** and keep the shipping function;
- planted logs recovered on every seed on the `GF(2^41)` and `GF(2^83)` test
  curves under the new walk, with measured expected iterations within the
  r-adding constant of the class-count prediction (`√(πn/2) / √262 × (1 + 1/(2H))`
  to first order);
- zero distinguished points failing verification, identical DP rate to the
  shipping walk within noise (the predicate is unchanged).

Inadmissible: changing `dpWeight`, changing the class structure (dropping
negation to save the `ε` computation costs `√2` in expected iterations, more
than the ~15% it would gain in rate), counting an ALU/clmad estimate as a
measurement, or quoting the rate without the correctness rows.

*Status (§6):* rows two and three met on the 4500 (300/300 device reports
re-walked per binary, 64/64 planted logs on `GF(2^41)`, DP rate within 0.12%);
row one is not decided: the 6000 could not be allocated, and the extrapolated
15.9 – 17.0 B/s straddles 17.0.

### 4.5 First measurement before any code

Re-run `benchmarks/clmad-price/probe.cu` on the 6000. Every clmad ceiling in
this note rests on the 4500's `1/37.7`; a 20% difference either way moves the
top of the achievable range by the same 20% and does not move the answer to
the 28 B/s question, but it decides whether the `R + ε·σᵏ(T_h)` walk is worth
17 or 20.

*Done on the 4500 only* (§3.1): `1/37.9` per `CLMAD.lo`, 63.9 slots per
64×64 product. The 6000 re-measurement is still owed.

## 5. What 28 B/s costs

Per GPU the floor is ~25 B/s at 100% of the carry-less unit and ~23 B/s at
90%. The route to 28 B/s of ECC2K-130 updates is therefore two RTX PRO 6000s at
≥ 14 B/s each — which the shipping kernel already does — or a `g7e.12xlarge`'s
four at the measured rate for 58 B/s. Per-GPU work above ~19 B/s buys nothing
the second GPU does not buy cheaper.

## 6. Measured

The walk of §4 is built behind `WALK_TABLE=1` (`include/tablewalk.h`,
`include/packedtablewalk.cuh`, host reference in `include/ref.h`,
`aws/campaign.json` field `walk`, which enters the campaign identity so the
two functions can never share a corpus). Everything below is from
[`benchmarks/table-walk/comparison.json`](benchmarks/table-walk/comparison.json)
and the raw logs beside it. Host: one **RTX PRO 4500 Blackwell Server Edition**
(82 SMs, 2.415 GHz maximum, **165 W limit**) on a `g7.2xlarge`, CUDA 13.3.1,
`sm_120`, the audited preset knobs. The 6000 was not available: every `g7e`
launch in `us-east-1` and `us-west-2` returned `VcpuLimitExceeded` against the
quota the running worker fleet occupies, so **no number in this section is a
6000 measurement**; §6.3 extrapolates and says so.

### 6.1 Two design rounds

| round | selection of `h, k, ε` | static ALU slots | 4500 B/s | vs shipping |
|---|---|---:|---:|---:|
| 1 | eight `L`-bit-plane masked popcounts (§4.3 as written) | 2,147 | 5.493 | +7.6% |
| 2 | byte LUT for the phase sum (17 × 256 B), nibble LUT for the pivot (33 × 16 B), inverse-log LUT; one `LDS.U8` per byte or nibble of `x` | **1,831** | **5.552** | **+8.7%** at fixed power, **+16.4%** per SM-clock |

Round 1 was a **relabelling** in part: it moved the walk's work from the
Frobenius network into 72 quarter-rate popcounts per update, and the count it
tracked (nominal ALU) fell 15% while the priced cost fell 7.6%. Round 2 moved
the same work onto the MIO pipe (shared-memory byte loads at 9.2 lanes /
SM-clock, otherwise idle in this kernel apart from the σ masks it replaced)
and is **engineering** against the same floor. Shared memory per block:
48,508 bytes (131 × 8 × 9 words of table + 4.9 KB of LUTs), two blocks of
256 threads resident per SM as before.

### 6.2 Correctness

| check | shipping | table (LUT) | table + ALU squaring |
|---|---|---|---|
| device selection primitives vs host reference, 4,096 points (`test-table-walk-cuda`) | n/a | phase 0, pivot 0, sign 0, tag 0, addend 0 mismatches; cycle rule fired on 1,536 | same binary |
| 300 device reports re-walked by the host reference, `dpWeight = 50`, 6 × 16 steps | 300 / 300, 0 dropped | 300 / 300, 0 dropped | 300 / 300, 0 dropped |
| distinguished points produced, identical iteration count | 331,611 | 331,227 (−0.12%) | 331,227 |
| planted logs on `GF(2^41)`, 16 instances × 4 seed salts, scalar reference engine | 64 / 64 | 64 / 64 | — |
| mean iterations to solve on `GF(2^41)`, `ℓ = 549,756,390,943`, 128 walks, `dpWeight = 13`, 64 samples each | 164,200 ± 9,900 | 157,800 ± 7,500 | — |

The class-count prediction for `GF(2^41)` is `√(πℓ/2) / √82 = 102,600`
iterations before the distinguished-point tail; 128 walks at `θ = 0.0138`
add ≈ 9,300 and the walk advances in launches of 128 × 64 = 8,192, so both
walks carry the same ≈ 1.5× harness overhead over the bare prediction (the
sample standard deviation is half the mean, as rho's is). What §4.4 asks is
the *ratio*: table over shipping is **0.96 ± 0.07**, consistent with the
r-adding constant `1 + 1/(2H) = 1.06` at `H = 8` and excluding a degradation
above 1.10 at two standard errors. Raw runs in
[`benchmarks/table-walk/raw/small-curve-gf2_41.txt`](benchmarks/table-walk/raw/small-curve-gf2_41.txt).

`GF(2^83)` planted logs, the second curve §4.4 asks for, need ≈ 2^37 steps of
the scalar reference and were not run; the device path is exercised on the
full curve by the 300-report re-walk instead.

### 6.3 Rate

Three alternating repetitions per binary, `--bench --steps 1024 --launches 32`,
SM clock and power sampled after each:

| binary | B/s, rep 1 / 2 / 3 | median | SM clock (MHz) | power (W) |
|---|---|---:|---|---|
| shipping | 5.151 / 5.107 / 5.069 | **5.107** | 1965 / 1950 / 1920 | 154 / 154 / 159 |
| table, LUT | 5.593 / 5.552 / 5.510 | **5.552** | 1830 / 1822 / 1807 | 154 / 166 / 153 |
| table, LUT + ALU squaring | 5.568 / 5.522 / 5.491 | 5.522 | 1830 / 1822 / 1987 | 154 / 152 / 152 |

Both kernels run the 4500 into its 165 W limit, and the table kernel draws
more per clock (its issue slots are fuller), so the card clocks it 6.6% lower.
The fixed-power gain is therefore **+8.7%**; per SM-clock the shipping walk
does 0.0319 updates (31.3 SM-clocks per update) and the table walk 0.0372
(26.9), a **+16.4%** gain. The per-clock figure is the one that transfers to
a part with power headroom; the 6000 has 3.2 W per SM against the 4500's
2.0.

Pipe utilisation at 26.9 SM-clocks per update: the static ALU figure of
1,831 slots is 28.6 SM-clocks at 64 lanes/SM-clock, so the dynamic count is at
most 1,722 and the ALU pipe is at its limit; clmad at ≈ 40 dynamic (the slot-0
arm's six are executed once per batch) is 23.5 SM-clocks, 87%. Moving the
squaring to the ALU (−5 clmad, +107 slots) measured −0.5%, which is what an
ALU-bound kernel does. **The shipping kernel was co-limited by both pipes; the
table kernel is ALU-limited with the carry-less unit at 87%**, so the next
slot to buy is ALU again, not clmad.

**Extrapolation to the RTX PRO 6000 (not a measurement).** The audited 14.64
B/s at 0.0319 updates per SM-clock implies the 6000 sustained 2.44 GHz in that
run. If the table kernel holds that clock, 188 × 2.44 GHz × 0.0372 =
**17.0 B/s**; if it loses 6.6% of clock as it did on the 4500, **15.9 B/s**.
The §4.4 target of 17.0 sits at the top of that range: **undecided.** The
measurement that settles it is the same three-way job
([`benchmarks/table-walk/gpujob.sh`](benchmarks/table-walk/gpujob.sh)) on a
`g7e`, and until it runs, the
shipping function stays the campaign default (`aws/campaign.json`:
`"walk": "sigma"`).

### 6.4 Classification

| change | class | evidence |
|---|---|---|
| quarter-rate pricing in `kernel_cost.py` | accounting | §3.1; no kernel changed |
| table walk, bit-plane selection | relabelling → engineering | nominal ALU −15%, priced −7.6%, measured +7.6% |
| table walk, LUT selection | **engineering** | priced −21%, measured +8.7% at fixed power / +16.4% per clock, correctness rows all green, ratio to floor 1.7 (bounded below by 1) |
| ALU squaring | engineering, did not pay | −0.5% |
| 28 B/s on one GPU | unchanged: **no** | floor rows of §2 |
