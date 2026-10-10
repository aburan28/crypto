# Larger norm-one field controls (2026-10-09)

The requested expansion is to larger curves in the existing norm-one
weak family, with both larger primes and larger odd extension degrees
confirmed by the user. Preserve the original complete-prime census request through
p=199; a larger-field control does not complete that census or classify
zero trace rows. The theorem covers every odd extension with q=1 mod 4.

| Requirement | Execution and evidence |
| --- | --- |
| Larger p in the p^6 family | Point-count controls at p=257,509,1009,2003,8191,65537 and larger characteristics, with fixed seeds and recorded caps. |
| Degrees 5 and 7 over F_(p^2) | Controls over p^10 and p^14 at p=13,257,1009,65537, with separate receipts. |
| Fields near 192 and 252 bits | Attempt p=nextprime(2^32) and nextprime(2^42), degree 6; retain completion or timeout. |
| Original every-prime census through 199 | Resume from p=59 with persistent data and an independently managed worker. |
| Curve and isogeny verification | Norm and fourth root; two explicit independent rational halves on the target; point counts on source and target; Hasse and mod-16 checks. |
| Exact curve identities | Save finite-field modulus and coefficient vectors, source and target model records, and ICV1 identifiers for completed counts. |
| Class reach | Individual positive controls only; complete zero labels require the retained exhaustive class method. |
| Publication | Update existing PR 1588 with dated report, raw receipts, diagram, PDF, and relevant validation. |

The earlier temporary queue and p59 worker are absent at this round's
start. Its last committed snapshot remains historical evidence. Restart
p59; do not infer a completed result from its former RUNNING status.

Use one GP process per larger-field cell, 120 seconds for an initial
one-curve point-count probe and 256 MiB initial PARI stack. Check full
4-torsion by explicit halves rather than group-order factorization.
Large CM discriminant factorization is a separate stage and remains
unknown when it has not completed. Record actual Q and log2(Q); this
field-size quantity does not specify a prime-subgroup size.

Pre-register the full-range census work requirement C_p=(p^4+3p^2+8)/12
using the retained theorem and identity IDC1h1bc37ec40bcb8619. For p=257,
509,1009 this is 363,555,713; 5,593,645,151; 86,374,331,401 calls.
These are exact algebraic counts, not measured runtimes.

The proved prime-degree orbit formula is
C_(p,n) = [N+A+B-3+2*e*(n-1)^2]/(4*n), for odd prime n,
N=(p^(2*n)-1)/(p^2-1), A=(p^n-1)/(p-1), B=(p^n+1)/(p+1),
and e=1 when p^2=1 mod n, otherwise 0. Exceptional orbits are
classified by the proof; eight cyclic-group and 45 stabilizer audits
pass. The cleared sum has record IDC1h90d58cc0e0c48fe3.

All 17 paired cardinality probes and all 49 geometric checks completed.
The temporary October 8 worker was interrupted. The replacement service
uses a verified internal-volume runtime, finite RunAtLoad=true /
KeepAlive=false launchd configuration, and persistent receipts under
/Users/adamburan/Library/Application Support/crypto-iso1/census-20261009-run3.
The first external-volume service failed an OS access check; its errors
are retained. The provisional restart-on-exit service was stopped before
completion and replaced by the finite service; its status is retained.
No external-volume access permission was changed.
