# Linear Sinz replay of the first n19 sparse S3 archive

## Question and immutable input

PR #795's first and only labeled n19 sparse S3 growth attempt exported a CNF but its independent verifier was censored at 600 seconds inside a redundant large-group witness sweep. The sweep ran after exact comparison with the ordered Sinz clause stream. This successor asks whether replacing that quadratic sweep with a proved linear structural check lets the **same immutable archive** pass the remaining independent point and path semantics. It does not rerun the exporter or change its representation.

The sole input is `../rotated_s3_sparse_n19_20260925/evidence/final/producer`, from Actions run `36526326477`, label event `32052162480`. The archive manifest, release receipt, producer result, old verifier and old freeze are hash-pinned in `FROZEN.json`. The original archive replay must remain `CENSORED`; a new semantic replay is a separate result. No SAT, PDP timing, full ECDLP cost, or n131 transfer claim follows from a pass.

## Structural theorem used by the new checker

For primary literals \(p_0,\ldots,p_{k-1}\) and fresh auxiliaries \(s_0,\ldots,s_{k-2}\), the ordered block is:

- \(p_0\lor\cdots\lor p_{k-1}\);
- \(\neg p_0\lor s_0\);
- for \(1\le j\le k-2\), \((\neg p_j\lor s_j)\), \((\neg s_{j-1}\lor s_j)\), \((\neg p_j\lor\neg s_{j-1})\);
- \(\neg p_{k-1}\lor\neg s_{k-2}\).

For \(k=1\), the block is simply \(p_0\). For \(k\ge2\), these are exactly \(3k-3\) clauses. The first clause requires a true primary. If \(p_i\) is true, its forward clause and monotonic auxiliary clauses force every \(s_j\) with \(j\ge i\) true. Its backward clause (or the final boundary clause when \(i=k-1\)) and monotonic clauses force every \(s_j\) with \(j<i\) false. Two true primaries \(p_i,p_h\) with \(i<h\) would therefore force \(s_{h-1}\) both true and false. For exactly one true \(p_i\), the unique assignment \(s_j=[i\le j]\) satisfies all clauses. Thus the block has exactly one satisfying auxiliary extension iff exactly one primary is true.

The checker verifies positive, distinct primary and auxiliary variables, then compares every archived clause against this independently generated sequence once. The existing schema check independently verifies global variable numbering and disjoint groups; the streaming DIMACS reader rejects extra, missing, malformed, or out-of-range clauses. Complete truth tables for \(k=1,\ldots,5\) and two corruptions of the actual S6 block (middle-clause sign flip and deletion) guard the implementation. This replaces \(O(k^2)\) witness scanning with \(O(k)\) clause matching while preserving the same exact-one semantics.

## Replay and decision

Run Python 3.12 `test_linear.py`, the original `ci_replay.py --evidence` on the first archive, then `verify.py --produced` on its producer directory. `verify.py` keeps the original 600-second child wall and 512-MiB RSS caps, parsed transition check, all 117,649 signed point tuples, 58,825 archived path comparison, 33 target checks, and five geometric mutation controls. Its only semantic change from the archived verifier is the one-hot block check. Source, input, and archive hashes must match `FROZEN.json` before replay.

A complete `PASS` requires every check and control. A wall or memory cap is `CENSORED`; a hash, clause, path, point, target, or control mismatch is `FAILED`. The replay's wall and memory values describe verification overhead and are excluded from any attack-cost comparison. On `PASS`, the next separately frozen gate may time same-target SAT solvers and independently lift every model; the n19 result alone does not establish feasibility at n131.
