# Proof map: cleared χ chord resultant

The readable derivation is in [RESULTS_V1.md](RESULTS_V1.md). The pinned
Lean 4.34.1 project is
[`inverse_velu/lean`](../../portable_order_trend_20261015/inverse_velu/lean/),
and the new module is
[`ChiChord.lean`](../../portable_order_trend_20261015/inverse_velu/lean/ChiChord.lean).

| report claim | Lean declaration | exact hypotheses | status |
| --- | --- | --- | --- |
| Cleared chord polynomial vanishes after substituting the line, χ-coordinate and source-curve equations | `ChiChord.eliminated_of_line_and_curve` | all integer `m,d,x₁,y₁,x₂,y₂,n,a,b`; `m(x₂−x₁)=y₂−y₁`, `n=dy₂`, `y₂²=x₂³+ax₂+b` | `lake build ChiChord` passed; axioms `propext`, `Quot.sound` |

The theorem's `a,b` specialize to `−21,37` in the run. It proves the
algebraic implication conditional on the line and point equations. It
does **not** define `F_(211²)` in Lean, prove that the specific Rust
line branches cover all exceptional cases, or prove root-finding
completeness and multiplicity over every field. The exhaustive 231-case
check and independent Rust checker cover the fixed public instance,
not those universal formal obligations. A faithful finite-field and
curve-level theorem remains open in the existing
[`problems/nonrational-isogeny-fiber-dlp-cost-gap.md`](../../../../problems/nonrational-isogeny-fiber-dlp-cost-gap.md).

The project pins `leanprover/lean4:v4.34.1` in `lean-toolchain` and its
Lake manifest. `#print axioms` reports only Lean's standard `propext`
and `Quot.sound`; the theorem has no `sorry` or custom axioms.
`lean4checker --fresh` was unavailable on this host, so that extra
independent kernel pass is not reported as done.
