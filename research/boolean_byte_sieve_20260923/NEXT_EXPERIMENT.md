# Prospective direct schedule construction, with a cost-bound gate

This is a new hypothesis, not an implemented or measured gain. The current
experiment remains rejected under its universal criteria. Its fresh holdouts
must become regression inputs in any successor; they cannot be tuning data for
another claim against the same holdouts.

The compiled full-word schedules are the strongest family in this round, but
the current constructor first creates a dense global syndrome representation
and then derives the schedule's low-coordinate images and high differences.
The possible degree of freedom is to build those required objects directly from
the original canonical monomials. This preserves the requested coefficient-aware
contract while potentially avoiding a redundant intermediate representation.

Before implementing a constructor optimization, measure construction and scanning
separately on discovery inputs, preserving outcomes, models and semantic work.
Record the instrumentation overhead against the retained complete-solve control.
The current `enumeration_ns` combines setup and scanning and cannot supply this
split. If the unchanged scan already costs T_scan and the strongest complete
reference costs T_ref on a fixture, even free construction cannot provide more
than T_ref/T_scan on that fixture under the stated unchanged-scan assumption.
Do not pursue a construction-only dramatic claim if that bound is insufficient.

If the cost split justifies the experiment, compile exact contributions as follows:

- Constants, low linear terms and low-low quadratic terms contribute to the
  initial 16- or 64-point full-word block.
- High linear terms contribute to their scheduled first differences.
- Low-high terms contribute to the corresponding fixed cross-coordinate image.
- High-high terms contribute to lower-difference updates and, when adjacent in
  the declared order, the required initial neighbor term.

Every contribution uses XOR. Repeated terms must cancel; zero and constant
generators, degree drops and changed coefficient vectors remain valid. No
fixture witness, target-specific precomputed answer, hash-only equality or
unpriced schedule is available to the solver.

Compare each directly compiled object against construction through the full
syndrome on exhaustive small inputs and dense word-boundary inputs. Keep the hot
scan shared where possible so a changed memory layout or compiler specialization
is not silently attributed to construction. If the hot representation changes,
that is a separate mechanism with a separately charged comparison.

Any accepted candidate needs the retained current methods, the same quiet controls,
fresh inputs, paired randomized/reverse method orders, complete cold costs and the
unchanged >2x criterion. Preserve individual exceptions, missing completions and
censored costs. Confirmation requires unchanged timed source on unused holdouts.
Production, full index-calculus, calibrated-operation and rho costs remain null
until measured. No asymptotic or cryptanalytic result is proposed here.
