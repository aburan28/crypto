# Long-XOR interface contract

This is a new thin wrapper linked to the frozen modified CMS functional-primary core. Only the XOR length check changes from3 to16,384, and the unsupported implied_by interface is disabled. Native configuration readback remains exact. No internal Gaussian configuration or core is changed.

There are8 primary bits,256 declared AND outputs g_j=b_(j mod8) AND b_((j/8) mod8), and declared output z. Each AND gate has its standard three exact clauses. The long parity equation contains all265 declared variables and has RHS0; z is forced0. Each ordered pair of input bits occurs four times among g_j, so their total parity vanishes. Thus accepted primary models are exactly the128 even-parity8-bit vectors. All declared auxiliaries are uniquely determined, and blocking only the complete8-bit primary vector enumerates each accepted model once. UNSAT or a contradictory blocking clause is required for complete enumeration. UNKNOWN/errors/deadlines do not certify completeness.

The independent oracle enumerates256 primary vectors by integer popcount; replay checks every declared AND coordinate and parity directly. A separate empty parity row with RHS1 must be rejected. This tests long-row handling and native internal cut-variable interaction only. It does not encode the actual affine-cone norm predicate, prove original root coverage, establish novelty, or measure acceleration.

Review: constants, empty contradiction, repeated input pairs, uniqueness and full-coordinate model replay are retained. No use of pd_implied is permitted. Time includes setup, solving, replay and serialization. The outer supervisor charges child startup and completion. Late scientific completion must be demoted by the supervisor.

Next: full compressed-norm tiny complete-set controls, then a separately frozen hard attempt with native matrix statistics. No unchanged retry of the old hard timeout.
