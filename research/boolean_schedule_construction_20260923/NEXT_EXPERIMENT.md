# Prospective necessary affine fibers from equation-space projection

No complete solver or speed ratio for this method is implemented. The structural
census uses only the 24 discovery fixtures, not prior or new holdout timings.

Choose a fixed low-coordinate block I of size k. For each quadratic monomial
involving two coordinates of I, collect its equation-coefficient vector q_ij.
Let W be their span. Construct a fixed linear projection P with kernel W. Then

    P F(y,z) = P b(z) + sum(i in I) P a_i(z) y_i

is affine in the low coordinates y for every outside assignment z. This removes
the independent-set restriction: interactions inside I are annihilated in equation
space rather than forbidden in the variable graph. All coefficient operations
remain exact XOR operations, with fixed projection construction charged once per
cold solve and every changing coefficient update charged during the scan.

The projection is only a necessary condition. An original solution satisfies it,
but a projected solution can fail the original system. For example, projecting
`x*y+1` along its sole quadratic coefficient gives the zero system, while only
x=y=1 solves the original equation. A valid solver must retain and check every
original equation. If the affine fiber has free variables, it must examine all
remaining candidate extensions before rejecting the outside assignment.

The exact census covers k=4,5,6 on n12/16/20/24. At n24, the low-quadratic span
has ranks 3–5, 6–8 and 9–12 respectively. The corresponding quotient dimensions
are 21–23, 18–20 and 14–17. These are dimensions of the equation-coordinate
quotient, not counts of independent conditional equations or predicted rejection
probabilities. No uniformity or independence assumption supplies a speed claim.

Required bounded experiment:

1. Freeze the low-block choice and projection algorithm on discovery inputs.
   Charge rank computation, coefficient projection, storage and all setup.
2. Compare every projected coefficient with direct polynomial restriction.
   Include cancellation, zero projections, rank drops and constant equations.
3. Solve the necessary affine subsystem exactly. Preserve inconsistent and
   rank-deficient fibers; retain an explicit basis for free-variable recovery.
4. Verify all surviving extensions on the original system. Caps remain UNKNOWN;
   projected consistency never substitutes for an original solution.
5. Compare complete cold cost against the current compiled full-word kernels,
   prior independent-set fibers, quiet controls and retained methods. Use fresh
   holdouts and the unchanged dramatic-gain criterion if discovery merits a run.

The constructor-only cost diagnostic does not rule out this mechanism: it changes
the scan's mathematical work. Conversely, this structural possibility does not
repair or supersede the current cost diagnostic's failures and inconclusive cells.
Production, calibrated-operation, full-IC and rho costs remain unmeasured.
