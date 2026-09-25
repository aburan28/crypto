# Prospective implicit consistency polynomials

This is an unimplemented mathematical hypothesis, not a performance result. It
does not change the frozen projected-fiber run or authorize reuse of its holdouts
for tuning. The fixed scan above recomputes affine consistency as its coefficients
change. A different representation might express some of that consistency as
Boolean polynomials in the outside variables.

For the necessary affine system A(z)y=b(z) with k low variables, choose any k+1
equation rows R. Every original solution satisfies

    g_R(z) = det([A_R(z) | b_R(z)]) = 0.

Each column of A is affine in z and b is quadratic, so g has ordinary degree at
most k+2. Reduction in the Boolean quotient can lower its degree or erase it;
duplicate products cancel by parity. This identity follows by evaluating the
polynomial determinant at each assignment. No independence assumption is needed.

The condition is only necessary. A rank-deficient A can make every chosen minor
zero even when the affine system is inconsistent. A zero polynomial is therefore
an ineffective filter, not a proof of solvability. The complete affine recovery
and original-equation checks remain required for every surviving assignment.

Before attempting a solver:

1. Use only the retained n12/n16 discovery fixtures. Fix k=4 and a deterministic
   row choice among non-pivot equation coordinates before inspecting its yields.
2. Construct a small capped set of minors using exact Boolean multiplication
   and XOR cancellation. Retain zero minors, degree drops, duplicate minors,
   compilation work, maximum support and retained storage. Exceeding a cap is
   an incomplete construction, not a small polynomial.
3. Compare the polynomial value with an independently evaluated determinant on
   every bounded outside assignment. Check necessity against direct original
   solutions and include a rank-deficient counterexample to sufficiency.
4. Measure full cold construction and evaluation against direct consistency
   checks on the same discovery inputs. Count overlapping rejections exactly;
   do not add individual filter rates or assume independence.
5. Only if that diagnostic justifies a complete candidate, freeze it and compare
   complete solves against the retained fastest references with new holdouts.
   Charge polynomial construction, all evaluations, affine recovery and original
   checks. Reusing coefficient-dependent output from another system is invalid.

This returns to the symbolic-product question with a distinct contract:
determinant coefficients depend on products of changing input coefficients, not
only on a fixed linear routing of one coefficient vector. A support envelope must
represent those products and cancellations explicitly. Neither a degree bound nor
a small symbolic support establishes a fast complete algorithm.
