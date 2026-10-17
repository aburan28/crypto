# Fresh prime-field elliptic-curve IC attempt

This is an actual bounded native C++ Semaev/Gröbner relation pipeline, not the
pinned C backend's multiplicative-group index calculus. It runs on the fresh
fixed P/Q fixtures without receiving target scalar witnesses.

For each generator: sample twelve independent unknown-log sign orbits on the
source curve with getrandom/rejection sampling and transport the same base
through each chosen map. Check every base point and exact prime order. Freeze
twenty nonzero uniform known setup coefficients t, paired across curves. Attempt
twenty R=[t]P relations per base, reduce valid signed two-sum coefficient rows
modulo r, record rank, and recover base logs only if rank reaches twelve.

For each actual unknown-scalar target Q, compute the two-sum Semaev ideal
(f3(X,Y,xQ), product(X-xFi), product(Y-xFi)) over F_p. Solve with native grevlex
Buchberger, sugar-prioritized pairs and product criterion. Retain counters,
critical-pair events and returned bases, extract roots from that basis on the
factor-base grid, then validate actual signed point sums. Run two positive
controls per base: F0+F1 and [2]F0. The independent signed-point lookup is built
only for verification and never passed to the polynomial solver or relation
collector; its setup cost is retained separately. Grid-variety equality and
Buchberger/input membership checks are outside the measured GB interval.

Record construction, GB wall/process CPU, root/relation and verification costs;
all native field-operation and polynomial degree counters; relation yield,
rank, failures and recovery status. The maximum actual polynomial degree and
maximum processed-pair sugar degree are implementation metrics. Define dreg
separately as the first zero Hilbert-function degree of the leading homogeneous
ideal (X^12,Y^12,X^2Y^2); this is not CM regularity of a homogenized ideal.

The twenty setup probes per base are an explicit initial attempt cap. Rank
below twelve gives insufficient_setup_rank, not a solved logarithm or a
mathematical impossibility. Each actual target's PDP/degree metrics are still
measured. The full IC recovery requirement remains open after a capped attempt
unless the candidate scalar is independently verified against Q=[k]P.

All source and IC inputs are hashed before the first probe. Source changes
create new attempts; incomplete records are preserved. CPU-only L0 with
concurrent native solver work is explicit; no GPU, isolation, complete-cost
speedup or external-fastest claim is made. Cross-process/external GB replay is
an additional audit obligation beyond the included exact grid/point checks.
