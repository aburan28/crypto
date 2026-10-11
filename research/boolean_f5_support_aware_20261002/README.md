# Support-aware exact Boolean F5 batch continuation

This is the frozen [protocol](PROTOCOL.md) and native development implementation
for the next bounded experiment
after the [rejected graded batch discovery](https://github.com/aburan28/crypto/pull/1236).
The parent result measured a real 1.78–1.84× complete-batch improvement on
public generated systems, below its universal 2× gate. Its fresh holdout
was not used. This protocol changes the *mathematical support contract*:
when the fixed high prefix is fully present but lower columns are absent,
the same row transform can be applied before projecting onto the exact
source column layout and resuming M4RI. The candidate must still return
the identical ordered Boolean F5 output and charge every fallback.

The native worker contains the same-binary rejected anchor, guarded support
projection, direct output unpack and a retained quadratic row layout with
scratch arena. It recomputes the Boolean F5 criterion on every assignment.
Development tests at n=12/16/20/24 compare selected rows and ordered output
exactly with the inherited F5 implementation, including sparse-support
cases the anchor sent to fresh F5. A local ignored timing test uses only
development seeds and is not a qualified performance result.

On a committed clean Linux x86-64 AVX2 checkout, use
`bash research/boolean_f5_support_aware_20261002/run.sh discovery NEW_DIRECTORY`.
The script reserves logical CPU 2's full SMT sibling core, retains
readiness and isolation receipts, charges complete cold batches, and seals
native replay without overwriting an attempt. The new discovery and
holdout seeds are fixed in `protocol.json` and have not been executed. No
verified speedup, full index-calculus cost or Pollard-rho comparison is
claimed here. There is no curve or key input.
