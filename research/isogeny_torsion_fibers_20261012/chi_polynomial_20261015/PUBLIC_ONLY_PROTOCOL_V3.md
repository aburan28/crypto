# Public-only same-target repeat: frozen before implementation

The V2 root-oracle and collector interfaces saw only public points, but
the V2 runner read `answer_scalar` from the fixture for the framework's
final assertion. V3 removes that scalar from the *entire solver process
input*. Retain V2 source, raw records and replay as historical evidence;
this is a new run, not a replacement.

## Frozen comparison

Use exactly the four V2 `F_(p²)` cells `(p,r)=(61,19),(271,61),
(787,211),(3319,811)`, the same `G,Q`, signed projected base, line,
cofactor, three round seeds, trial cap `min(20r,10^6)`, 10⁷ direct-root
polynomial-call cap, 512 MiB shared-process RSS cap, and all 11 arms in
the same rotated order. Freeze a public JSONL by removing only the
`answer_scalar` field from every V2 fixture object; the run spec pins
both old and public content hashes. The public-only executable must
reject an input object containing `answer_scalar` and must not open the
private fixture or answer file. The independent replay process alone
may read the private fixture.

Add a native framework entry point that retains every factor-base,
oracle, relation, matrix and source-point verification phase but omits
the optional planted-scalar equality assertion. Its `verified` report
flag must still require `[k]G=Q`, already checked in the relation
collector. Keep the planted entry point's behavior unchanged and test
that both entry points agree on a known public synthetic instance.

First check that public and private fixtures differ *only* by the four
removed scalar fields. Run the public-only solver once on the frozen
four-cell, three-round workload; preserve source/binary hashes, command,
stdout, stderr, timings, primitive counts and failures in new `*_V3`
files. The independent checker must check every 11-arm slot, verify
`[k]G=Q` with separate arithmetic and optionally compare to its private
answer fixture. Any missing slot, changed input, scalar mismatch or
failed arm outside the inherited raw-pair cap exhaustions censors the
affected result. Stop after one run and one replay; do not rerun for a
favorable timing. Compare V3 to V2 as a custody check, not a speed
estimate. The method reference remains same-target direct negation rho;
portable counted cost, isolation and independent-host performance replay
remain separate gates.
