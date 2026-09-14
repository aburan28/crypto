# What the block experiment proves, and what to change next

The soundness question is settled for this relaxation: a saved dual separator
proves an entire z block impossible at fixed b. The broader cost hypothesis
is not established. The new bilinear representation improves on the old hybrid,
but indiscriminate pruning regresses at n18. At n30,d8 it reduces work against
the unpruned circuit, yet the direct S3 table remains the stronger enumeration
baseline, especially when setup is shared across eight targets.

## A new constructor-specific floor

The frozen `Circuit` constructor uses exactly

    3 + 8*d + 3*d²

field multiplications per b: three outside the coordinate loop, eight for each
pair of linear columns and mixed-row scalars, and three per mixed entry.
Both controls construct this circuit before attempting to reject the root
block. The expression is instrumentally checked by analyze_batch.py, and the
source hashes freeze the implementation to which it applies.

The complete n30,d8 batch constructs 2,048 such circuits. Its constructor
alone therefore needs 2,048*259 = **530,432 multiplications**, before rank tests,
inversions, extraction or verification. The matched S3 batch uses **207,401
multiplications in total**. This is a new component lower bound for the frozen
constructor, not the predecessor's seven-multiplications-per-leaf bound and
not a universal algorithmic lower bound.

Even free, perfect block pruning leaves the constructor 2.5575 times above the
S3 multiplication total. To be 20% below that S3 component, more than 68.71%
of these constructor multiplications must disappear before paying any other
cost. Exact operands and the computed ratio are in batch_summary.json. This
does not calibrate additions, squarings, binary work or full-DLP cost.

## Cheap symmetry is valid but insufficient

The displayed bilinear formula gives Lz=Lw and M(z,w)=M(w,z). Copying the
identical linear columns saves 2*d multiplications per b. Computing only the
upper triangle of M saves 3*d*(d-1)/2. At d8 those two changes reduce the
constructor to 159 multiplications per b, still **325,632 per batch**, above
the entire measured S3 multiplication count. The formula is an algebraic cost
projection, not a timed implementation result. Symmetry alone cannot cross
this component boundary with the remaining architecture frozen.

## A more substantive next hypothesis

Factor the mixed map using basis-only products. For basis elements e_i,e_j,
write

    T_ij = e_i*e_j,
    U_ij = T_ij²,
    W_ij = e_i*e_j² + e_i²*e_j,
    Q_ij = e_i²*e_j⁴ + e_i⁴*e_j².

Then the exact mixed entry is

    M_ij = Q_ij + (r²+h²)*U_ij + b²*r*W_ij + b²*r*h*T_ij.

Substitute h=b²+b+r to obtain

    M_ij = Q_ij + b⁴*(U_ij+r*T_ij)
                   + b²*(U_ij+r*W_ij+r²*T_ij) + b³*r*T_ij.

The first table is independent of target and b. For a fixed target, the
remaining fixed-scalar multiplications are binary linear maps on b⁴,b²,b³.
This suggests building reusable maps and updating the circuit columns as b
changes, possibly reusing elimination data as well. It moves work to table
construction, memory and XORs; none of those costs may be hidden. It is an
algebraic reorganization, not yet an implemented saving or a new exponent.

The next experiment should test this coefficient reuse against the frozen
unpruned and pruned circuits and the direct S3 table, charging table setup
both cold and over a named batch. Keep n18 regressions visible and include
fresh n30,d9 inputs: all complete enumerations timed out there in the current
three-second panel. Reject the candidate if field savings are offset by
elimination, binary work or memory. Increasing the timeout alone would not
resolve the cost boundary. The original three-size, calibrated 20% goal
remains the acceptance criterion.
