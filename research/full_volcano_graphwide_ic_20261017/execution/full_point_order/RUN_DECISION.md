# Before primary execution

Question: can the existing k=1 finite-field diagnostic retain full prime-power
images and recover the complete verified order of each selected public C point,
with matched direct PH controls, across every retained prefix model?

Success means every requested row is retained, maps are checked for exact image
order, answers are verified in the original cyclic subgroup, all actual work and
relation yields are recorded, and non-clock fields reproduce. A failed image
search stays a bounded failure; it does not prove mathematical nonexistence.
The unreduced sqrt(M) generic floor is not the relevant comparison for these
smooth orders. No performance win is a success condition or an admitted claim.
The direct reference is the explicitly implemented EC PH/BSGS control, not a
fastest-external or independently implemented solver.

Five native validation tests pass, including exhaustive scalar controls at
p=109/113/127, destination generator/digit controls, original C targets, wrong
orders and a rank-two outside-subgroup target. The first compilation's incorrect
counted-negation call was fixed before tests; that failure log is retained.
No original order-7 implementation or frozen survey binary was changed.
