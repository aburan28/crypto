# Stage 188 finalizer correction

The first charged final-verifier dry run used a draft result under
`development/`.  The verifier treated the draft's immediate parent as the
stage root and therefore tried to scan `development/development`, returning
exit 2 before checking the result.

The failed receipt is retained at `development/verification-dry-run/`.  The
additive correction discovers the nearest ancestor containing both
`PROTOCOL.md` and `development/`, so the same verifier supports a charged draft
below `development/` and the final result at the stage root.  No backend run,
raw metric, ratio, threshold, decision, or earlier verification is changed.
