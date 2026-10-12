# Prospective exact linear-fiber solver

No solver or cost claim is implemented by this note. The structural census uses
only the twelve discovery inputs from resource_probe_02. Primary holdout timings
must not determine a successor's thresholds, ordering or selection rule.

## Observed mathematical structure

Construct the union interaction graph: vertices are variables and an edge joins
i and j when x_i*x_j has nonzero coefficient in at least one original equation.
Parity is resolved within each equation before building the graph. An independent
set I has no quadratic term involving two of its variables. Thus, for every fixed
assignment z of the other variables, all original equations are affine in x_I:

    A(z) x_I = b(z).

The exact discovery census finds maximum sizes 2–3 at n16 and 4–5 at n24. Its
meet-in-the-middle algorithm uses a right-subset dynamic program and enumerates
all independent left subsets. Both sides have at most twelve variables; the
recorded optimum is finite and exhaustive, not a greedy lower bound. Source,
input-manifest binding, graph, explicit selected variables and work counts are
retained in linear_fiber_census_01. Python tests compare every graph through five
vertices with direct full enumeration, then replay every recorded case.

This establishes the existence of linear fibers containing 4–32 Boolean points
on those inputs. It does not establish a faster solve. In particular, the current
SIMD method already evaluates sixteen points together and can be cheaper than
building and reducing a small linear system.

## Required implementation and comparison

1. Freeze a deterministic variable-selection rule, and charge graph construction
   and exact selection inside every cold solve. Any simpler heuristic needs its
   own disclosed rule and cannot claim maximum independent-set size.
2. For each outside assignment, update A(z) and b(z) exactly. A(z) is affine in z;
   b(z) is quadratic. A Gray update is possible, but every coefficient change and
   setup cost belongs in the total.
3. Decide membership of b(z) in the column span of A(z), recovering a valid x_I
   when consistent. Rank can change, columns can vanish, and the zero matrix can
   be consistent or inconsistent. No full-rank or stable-pivot assumption is
   permitted. A missing scheduled pivot requires exact repair or full reduction.
4. Check recovered assignments on all original equations. Only exhausting every
   outside assignment permits UNSAT. A cap before exhaustion returns UNKNOWN.
   Distinguish enumerated outside assignments, linear systems solved and possible
   extensions proved impossible; do not label all of them as evaluated points.
5. Compare the coefficient updates and solution sets with direct restriction and
   exhaustive low-variable enumeration on small systems. Include rank-deficient,
   multiple-solution, constant, cancellation and zero-column cases.
6. Freeze a new matched protocol containing all current methods, including direct
   SIMD, both SIMD hybrids, block transport and initial restriction. Keep prior
   failures as regressions and use unused holdouts. Charge every phase and keep
   calibrated-operation, production, full-IC and rho costs null when unmeasured.

If |I|=k, there are 2^(n-k) outside assignments, but each carries a changing linear
system. The prospective ratio depends on the complete cost, not on the factor
2^k alone. This is a concrete degree of freedom beyond extracting the original
affine row span; neither this note nor the structural census establishes a gain.
