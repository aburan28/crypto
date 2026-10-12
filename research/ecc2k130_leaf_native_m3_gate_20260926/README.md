# Conditional native-leaf m≥3 gate

[`PROTOCOL.md`](PROTOCOL.md) isolates the unmeasured m≥3 native-leaf versus
exact source-pullback comparison and its complete cost ledger. The protocol
is held until an implicit full-point producer passes semantic and capacity
gates. It includes no new result or speed claim.

The read-only admission helper computes a **necessary upper bound**, for
example reproducing the 32-target, 16-point m=2 ECC2K-130 smoke ceiling of
[#753](../ecc2k130_factor_base_replication_20260925/README.md):

```sh
python3 research/ecc2k130_leaf_native_m3_gate_20260926/admission.py \
  --unordered-size 16 --arity 2 \
  --target-space-size 680564733841876926932320129493409985129 \
  --targets 32 --required-rank 1
```

This gives `support_cardinality_upper=136` and
`expected_supported_probes_upper_fraction=4352/680564733841876926932320129493409985129`.
An m≥3 proposal must supply its own exact physical point counts, uniform
target space, total probe cap and independent rank requirement. Ordered
rotated slots use `--ordered-slot-sizes`; their tuple product must not be
replaced by the unordered common-base formula.
