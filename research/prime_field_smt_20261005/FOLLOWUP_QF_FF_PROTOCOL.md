# Follow-up: native finite-field SMT encoding

Status: preregistered after the frozen `QF_BV` run returned `unknown` at
120,000 ms and before any `QF_FF` solver run.

Compatibility addendum, recorded before the frozen `QF_FF` run: the first
toy control with the non-GPL static artifact exited immediately with
`cvc5 can't solve field problems since it was not configured with --cocoa`.
That failed control is retained. Native finite fields require CoCoA, so the
remaining controls and frozen run use the same cvc5 1.4.1 release and git
revision from the official `cvc5-Linux-x86_64-static-gpl.zip` artifact,
SHA-256
`d0b54324ec2129697975da8753767fd255309947b87832917d32a78da7d16666`.
The extracted executable hash is recorded in every receipt.

The first protocol and its input, target, factor-base bound, cvc5 release,
time cap, verification rules, and claim exclusions remain unchanged. This
follow-up changes only the arithmetic representation.

## Hypothesis

cvc5's native prime finite-field theory will avoid the bit-vector division
circuits created by `bvurem` and return a model for the frozen 20-bit case
within the same 120,000 ms cap. The model must pass the same independent
big-integer curve and group-law replay.

## Frozen encoding

- logic: `QF_FF`;
- cvc5 option: `--ff-solver=split`, selected before any run because the
  cvc5 1.4.1 theory reference specifically identifies it as the solver for
  field equations that encode bit decomposition;
- field sort: `(_ FiniteField 630043)` for the frozen case;
- `x1`, `y1`, `x2`, `y2`, and `lambda` are field elements;
- each abscissa is linked to a little-endian vector of field variables
  constrained to `{0,1}`;
- the weighted bit sum equals the abscissa;
- a lexicographic Boolean formula constrains the represented integer below
  the unchanged bound `4096`;
- ordinary addition and doubling use the same denominator-free equations as
  the `QF_BV` query;
- ordinary addition uses an exact bitwise `x1 < x2` predicate, retaining the
  same symmetry breaking.

Because every represented integer is below `B <= p`, the bit-to-field map is
injective. The bound and ordering predicates therefore have their ordinary
integer meanings rather than meanings modulo `p`.

## Runs and stopping rule

Run the planted SAT and complete-UNSAT toy controls first, then run the
frozen cross-repository case exactly once with 120,000 ms. A timeout,
`unknown`, malformed model, or failed independent replay is retained and
falsifies the hypothesis. No alternate bit ordering, solver option, bound,
or target is tried after observing the frozen result.

This remains a decomposition-stage diagnostic. Pollard-rho ratio,
end-to-end `S`, and speedup remain `null`; no scoreboard row changes.
