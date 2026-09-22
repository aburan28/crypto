# Round 0022 — what round 0021 got wrong, and what the ladder says instead

Round 0021 reported the first measured crossing of rho by this collector:
`n23a1` at 0.831 and `n37a0` at 1.533, both bands clear of one. The first
number survives this round unchanged. **The second does not, and neither does
the reasoning that produced it.**

This round exists because two checks round 0021 did not run both fail.

## 1. rho was charged while it was still running

`round21_crossover.py` verifies that the IC arm completed, and then measures
both arms:

```python
report = run(worker, dict(job, mode='ic'), env)
if report.get('status') != 'complete':
    continue
...
ic_ir.append(instructions(worker, dict(job, mode='ic'), env))
rho_ir.append(instructions(worker, dict(job, mode='rho'), env))
```

`instructions()` parses callgrind's `Collected:` line. It never sees the
worker's JSON, so **it cannot tell a solved rho from one that ran out of
budget.** Re-running round 0021's own draw:

| cell | fixtures | IC complete | rho complete |
|:--|--:|--:|--:|
| `n23a1` | 64/64 | 64 | 64 |
| `n37a0` | 64/64 | 64 | **60** |

Four of the sixty-four rho runs at `n37a0` were cut off and charged anyway.
That is the campaign's own rule — an arm that does not complete has no cost,
and a budget exhaustion is never evidence about an algorithm. Round 0021
enforced it by hand at `n41a0` and missed it at `n37a0`, where the claim "both
arms complete at both cells" was asserted rather than checked.

### The direction of the bias, which I also got wrong

My first correction said the error inflated the ratio and that 1.533 was an
upper bound. **That was reasoning, not measurement, and it is false.** It
assumed a truncated rho is charged *less* than a finished one. The opposite
holds: a rho that returns `incomplete` has exhausted its restarts, so it did a
great deal of work and produced nothing, and its instruction count is *larger*
than a run that found its collision early. Charging it over-charged the
denominator and pushed the ratio **down**.

Measured, with every rho run required to complete: `n37a0` reads **1.763
[1.547, 2.009]**, against 1.533 published. **1.533 was a lower bound.**

## 2. The trial cap is protocol, and it is worth 14 ppm

Letting rho finish means raising `max_trials` from the frozen 4096 — the worker
itself accepts 65,536. The convenient argument for doing that freely is that a
cap which is never reached cannot be observed, so the raise would be additive
and every published measurement would read unchanged.

**That was predicted and it is false too.** `round22_budget_effect.py` finds
identical recovered logarithms and identical factor bases at every cell, and a
different instruction count at every cell on both arms, because

```rust
max_iterations_per_restart: cfg.max_trials,
```

makes `max_trials` rho's iterations-per-**restart** — how often it abandons a
walk and starts over, a parameter of the algorithm being measured — while the
candidate takes it at `TinyIc::new`.

So the exact-equality gate was right to reject it. What the gate could not say
is how much it matters, and that is the number that decides whether a ladder
measured at one cap can be read against a panel measured at another:

| cell | mode | Ir @4096 | Ir @65536 | Δ | ppm |
|:--|:--|--:|--:|--:|--:|
| `n19a0` | ic | 995,233 | 995,229 | −4 | −4.0 |
| `n19a0` | rho | 2,727,298 | 2,727,264 | −34 | −12.5 |
| `n23a1` | ic | 2,598,899 | 2,598,922 | +23 | +8.9 |
| `n23a1` | rho | 5,268,988 | 5,269,011 | +23 | +4.4 |
| `n37a0` | ic | 50,425,517 | 50,425,521 | +4 | +0.1 |
| `n37a0` | rho | 24,665,546 | 24,665,216 | −330 | −13.4 |

**|Δ| ≤ 14 ppm, and both signs.** A systematic protocol effect would push one
way, so the mechanism is worth separating: at `n23a1`, `max_trials` = 4096 and
9999 — same digit count, both far above anything reached — give the candidate
**bit-identical** instruction counts (3,277,995), while 65,536 differs by 19.
That share is the cost of parsing a longer JSON literal. rho does depend on the
value genuinely (4096 and 9999 differ by 46 out of 3.22M), which is the restart
policy, at 14 ppm.

`n23a1` is therefore carried at both caps as the ladder's anchor, and reads
**0.831 at both**, identical to three decimals. The ladder and the published
panel are comparable, with a measurement behind that rather than an argument.

## 3. The factor base was tuned toward the answer — by leaving it alone

Round 0021 wrote:

> **The factor base was not tuned toward the answer.** ... Bigger is
> monotonically worse ... The cheapest configuration is the one the sampler
> builds by default, and that is the one the headline uses.

The monotonicity is real on the eight-cell panel, where round 0019 measured it.
**It does not hold at either new cell**, and carrying it there was the round's
largest error, because the pair-table model puts the balanced optimum at
`F = (c·#E·t/k)^{1/3}` — which *grows with `r`*. A base pinned at the sampler's
default single batch while the optimum walks away from it measures the
candidate with a handicap that widens as the ladder climbs.

Sixteen fixtures a configuration, every report certified:

| orbits | `n37a0` IC/rho | `n43a1` IC/rho |
|--:|--:|--:|
| 8 (default) | 1.534 | 5.260 |
| 24 | — | 1.935 |
| 32 | **1.331** | 1.862 |
| 40 | 1.497 | — |
| 48 | — | **1.661** |
| 56 | 1.939 | — |
| 88 | — | 2.636 |
| 104 | 3.134 | — |
| 168 | — | 4.382 |
| 200 | 6.351 | — |

Both minima are interior. At `n43a1` the default costs the candidate a factor
of **3.2×**, on non-overlapping bands — 5.260 [3.669, 7.542] against 1.661
[1.328, 2.078].

Note which way this error ran: **reporting the default OVERSTATED the
candidate's loss.** It is the one direction a result that already goes against
the candidate will not be challenged on, which is exactly why the check had to
be run rather than inherited.
