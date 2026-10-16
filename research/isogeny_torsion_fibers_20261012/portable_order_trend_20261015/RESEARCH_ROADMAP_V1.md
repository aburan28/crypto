# Where the fiber and pairing experiments fit in the ECDLP program

The target is a **cold, single unknown scalar** on a public subgroup: return
`k` from `(E,G,Q=[k]G)` and independently check `[k]G=Q`. A candidate changes
the frontier only if its complete portable `S = work/sqrt(r)` and admitted
isolated wall cost beat the best same-target generic control over a credible
size trend. Multi-target reuse is a secondary, explicitly charged setting.
The [ecbench accounting contract](../../../docs/ecbench/README.md) and
[open cost problem](../../../problems/nonrational-isogeny-fiber-dlp-cost-gap.md)
remain the acceptance boundary. [Shoup's generic-group lower bound](https://www.shoup.net/papers/dlbounds1.pdf)
explains why an advantage must use structure or extra access; it does not
rule out a structure-specific algorithm.

```mermaid
flowchart LR
  A[Explicit curve and subgroup] --> B[Verified isogeny or pairing map]
  B --> C[New decomposition or finite-field solver advantage]
  C --> D[Relations or field log and scalar recovery]
  D --> E[Independent same-target replay]
  E --> F[Complete cold S and admitted wall vs rho]
  B --> G[Public equivalent-route search]
  G --> F
```

| Route | Current observation | Missing discriminating gate | What a positive gate would establish |
| --- | --- | --- | --- |
| Degree-three isogeny transport | Different-`j`, non-isomorphic source and destination pairs and exact three-point rational fibers passed independent finite checks. The [four-order panel](RESULTS_V2.md) and [selected order-1,201 addendum](order1000/RESULTS.md) solve every target. At order 1,201, full-image `S_group=1532.302` versus rho `4.071` at 64 targets; full-image construction enumerates the source. | Construct a useful destination factor base by a direct `χ`-coordinate resultant/root oracle without the source census. Exhaustively match its decompositions to an independent pair enumerator, then rerun complete recovery on the same public targets and controls. | A cheaper executable decomposition mechanism on this family. It would still need a complete cost crossover and broader orders. |
| Pairing to a multiplicative group | The [exact pairing transfer](../pairing_transfer_20261015/RESULTS.md) recovered both frozen scalars; torsion-input search dominated its two cold runs. | Price a reusable `H` and pairing base for a declared target count, together with an appropriate finite-field DLP, then compare with same-target direct rho. Keep cold and amortized rows separate. | Whether input-search reuse or a finite-field solver makes the transfer competitive for that field and target count. |
| Torsion-fiber algebra | `χ` and quotient-coordinate tests found structurally distinct support, but the first full-image support increase at order 229 did not persist in the three larger cells. | Prove and implement the direct polynomial oracle, including exceptional fibers and multiplicities, with a faithful curve-level Lean bridge where feasible. | Whether the coordinate geometry yields more useful relations per *fully priced* unit of work. |
| Hidden-witness isogeny route | No secret map or asymmetric construction is supplied by these public experiments. | Specify the withheld witness, designer setup and evaluation cost, public reconstruction cost, and cheaper equivalent routes, each on an exact family and subgroup; then price transfer and destination DLP. | A scoped assessment of asymmetric access, separate from a public universal attack. |

The most useful next algorithm experiment is the direct `χ` resultant. It
tests the exact bottleneck observed in the [full-fiber audit](reports/run-20261015T2335Z-fibers/REPORT.md):
building the full rational-image base by scanning `#E(F_(p²))` points in
every cold child. The selected order-1,201 audit counted 1,531,275 source
points, 510,425 rational images and 1,442,401 subgroup pair identities;
its full-image/destination-direct support was 122/128 keys. The oracle
must output a *complete, independently checked*
decomposition set, not just a faster lookup on planted relations. Its
reference is direct source BSGS and rho on the same 64 frozen targets,
with the existing seven-arm panel retained. [Semaev's summation-polynomial
construction](https://eprint.iacr.org/2004/031.pdf) supplies a public
algebraic control; the missing step is a useful `χ` specialization with
an executable root process, measured cost, and recovery.
[Joux and Vitse's cover/decomposition work](https://eprint.iacr.org/2011/020)
is another structural comparison, especially for composite-degree
extensions; it does not supply the missing `χ` oracle for these exact
quadratic-extension cells.

The selected fifth order extends the measured range from 229 to 1,201
(5.24×), with three cells clustered near 430 and varying field primes
and cofactors. This is a larger-order count observation, not a stable
exponent estimate. The next size experiment should keep the construction
and factor-base shape fixed where possible while varying subgroup order
more broadly. Preregister at least four geometrically spaced orders,
one unseen public target per cell plus a 64-target reuse panel, a common
512 MiB memory ceiling, and the same rho/BSGS controls. A cap or memory
failure is a censored cost result, not an absent curve or failed theorem.
This experiment is distinct from the completed five-order descriptive
count panel.

The [frozen price-sensitivity result](sensitivity/RESULTS.md) shows that
at 64 targets the omitted rho work would need to add 112.973–334.844
group-addition equivalents per observed step across the four earlier
orders to tie their full-image group-count lower bounds *if omitted
full-image work were free*. This is an algebraic threshold, not a
measured native price. The **pricing gate** is a pinned conversion from every counted native
category (field arithmetic outside group calls, square roots, hash/table
work and modular rows) to the same group-addition equivalent. The current
`S_group` and separate Fp axes in [TREND_V1.json](TREND_V1.json) are
portable, but `complete_portable_S` is explicitly unpriced. Instrument
missing operations, derive and freeze conversion intervals on at least
the measurement host and an independent host, and report a sensitivity
range. A componentwise loss across the priced axes may settle this
implementation's ranking without a fragile point estimate, while an
axis trade needs the conversion.

The **performance gate** is an interleaved L2-or-better session on a
reserved physical core, including rho A/A and exact resource receipts.
The [access check](ISOLATION_ACCESS_V1.md) identifies the current host
placement blocker. The [independent Linux replay](INDEPENDENT_HOST_PROTOCOL_V1.md)
checks portability and cross-host performance but cannot substitute for
that physical isolation. No canonical frontier or security claim follows
from the present shared-host timings.
