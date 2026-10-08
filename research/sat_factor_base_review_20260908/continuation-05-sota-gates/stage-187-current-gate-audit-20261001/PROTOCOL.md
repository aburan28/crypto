# Stage 187 protocol: current gate audit through dense-pair selection

Compose an additive seven-gate audit over immutable Stages 174–186 and the
current external reproduction thread.

The audit must:

1. Rehash every referenced `result.json` and `verification.json`.
2. Recompute the Stage 175–186 measured cost increment from unique components,
   excluding Stage 180's inherited Stage 179 five-column build and using the
   explicitly incremental fields in Stages 182, 184, and 186.
3. Add that increment to Stage 174's measured campaign lower bound. Keep
   complete campaign cost `null` because interactive compile/test work outside
   process meters still exists.
4. Record the final selected implementation: current five-column BlockTables,
   dense exact pair selection by default, quadratic selector at
   `F4_F2_DENSE_PAIR_SELECT=0`, and full M4RI retained only at
   `F4_F2_FULL_M4RI=1`.
5. Preserve every negative and accepted decision from Stages 175–186.
6. Snapshot the public external reproduction issue and count comments from
   users other than `aburan28`; no self-comment counts as external review.
7. Carry forward the seven gate labels without promotion unless source evidence
   changed. Licensed Magma, full IC versus automorphism-rho crossover, and
   unaffiliated reproduction/novelty review remain required.

The audit must state that dense pair selection is an internal single-target
engineering improvement, not Koblitz index-calculus SOTA.
