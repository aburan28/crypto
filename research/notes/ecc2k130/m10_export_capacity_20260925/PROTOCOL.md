# ECC2K-130 n131 m10 complete-chain representation-capacity preflight

Status: **frozen implementation proposal; no n131 complete-chain outcome**. This draft is stacked on the held #804 DIMACS gate, which depends on #802. The unequal input is copied byte-for-byte from held #784; balanced #778 is merged. `FROZEN.json` has `release_main_head: null`. The measured producer must not run until #802, #804, and #784 are merged, this branch is rebased on a fixed main head, the hashes and caps are frozen again, CI passes on that exact head, and the exact-head implementation receives independent review. A failed prerequisite or cap produces a STOP/censored receipt, never a result inferred from silence.

## Question and scope

Can the exact, complete ten-factor K0 addition-chain Boolean DAG, with compact implicit rotated factor domains, fit the frozen raw-DAG and deterministic DIMACS-size budgets? This measures **representation capacity only**. It does not run a SAT/ANF solver, sample a PDP target, estimate PDP yield, compute signed-column mapping/rank, or make an ECDLP/DLP claim. The generic and canonical-infinity target unit variants are counted solely to expose deterministic CNF length. No n131 CNF file is written.

The prior balanced screen [#778](https://github.com/aburan28/crypto/pull/778) measured n131, m10, ten 13-dimensional rotated slots: 7,977 physical F0 points, 3,988 nonzero signed columns, and necessary tuple count divided by q of 1.5329446704. The unequal screen [#784](https://github.com/aburan28/crypto/pull/784) measured dimensions `[14,13,13,13,13,13,13,13,13,13]`: 16,125 points in the high slot, 7,977 in each low slot, 8,062 normalized global nonzero signed columns, and necessary tuple count divided by q of 3.098750509. Both pass those *necessary* count/column gates, not a PDP-yield test. Their previously reported raw affine-chain floors of 1,178/1,179 variables are a different representation and are not substituted for the full-point circuit below. The small one-hot 7-point choice in [#811](https://github.com/aburan28/crypto/pull/811) is likewise not extrapolated to either n131 domain.

## Exact representation

Use the source curve `y^2 + xy = x^3 + 1` over `F_(2^131)`, the bit-polynomial modulus and normal element beta=3 in `INPUT.json`, and the already checked complete affine/infinity relation in [#802](https://github.com/aburan28/crypto/pull/802). For factor slot i and its frozen dimension d_i, create d_i Boolean selectors and wire the coordinate bits of

    x_i = XOR over j<d_i of (a_ij * beta^(2^(10j+i))).

The infinity flag is the constant false; its y coordinate is a free 131-bit word. The parent full-point relation enforces the on-curve equation on every operand. Thus each nonzero liftable x admits both rational sign lifts, a nonliftable x admits none, and x=0 admits the unique point (0,1). This is exactly the physical affine factor set for each slot, without one-hot enumeration or hidden point-choice table. The independent bit-polynomial replay checks beta's 131 conjugates have rank 131, balanced slots have combined coordinate rank 130 (conjugate 130 omitted), unequal slots rank 131, every individual slot has full rank, and six fixed masks in each slot reconstruct the same 131-bit x words.

Create eight full-point prefix roles S2 through S9, one full-point SUM role, and nine 131-bit existential slope words. Copy the complete relation nine times for F0+F1=S2, S2+F2=S3, ..., S9+F9=SUM, binding each shared point role to the same DAG nodes. The full-point roles each have a canonical infinity flag plus 131 x and 131 y bits. The primary input count must be 4,986 for balanced and 4,987 for unequal: `sum(d_i) + 10*131 + 9*(3*131+1)`. Nodes are hash-consed through #802's exact XOR/AND DAG; there is no partial-chain or branch omission. Under the node cap, #804's one-based Tseitin templates determine the exact generic and target-O DIMACS clause and ASCII byte counts. The target-O variant adds 263 unit clauses and is *not* submitted to a solver.

## Frozen gate and archived evidence

Per arm, independently: external wall cap 600 s, address-space/process-tree RSS cap 2 GiB, raw DAG cap 2,000,000 nodes (including constants and inputs), target-O ASCII DIMACS cap 256 MiB. The external parent watchdog comes from #804's `bounded.run_child`: it starts a separate process group, installs RLIMIT_AS before execution, polls process-group RSS, kills on the wall/RSS cap, and archives stdout/stderr and any child checkpoint. The child also sets RLIMIT_AS and an alarm as a backup; failure to install the memory cap is STOP. It fsyncs a checkpoint after primary wires and each completed edge, recording stage, counts, and prefix SHA-256. If the 2,000,001st node appears, `CENSORED_DAG_NODE_CAP` records stage and exact partial prefix; if a complete DAG is built but target-O bytes exceed 256 MiB, `CENSORED_CNF_BYTE_CAP` records the exact counted size. Neither censored status implies anything about PDP hardness, solver performance, or full uncapped DAG/CNF size. External time/memory stops preserve the raw and partial receipts and do not become negative feasibility claims beyond those caps.

The producer archives each arm's hashes, counts, exact byte sizes if complete, resource measurements, progress, and a manifest. A separate verifier process rebuilds the same DAG under the same cap and checks counts/prefix/checkpoints, while its independently written clause-length arithmetic replays both DIMACS byte counts. It independently regenerates the normal basis and fixed-mask coordinates. Because the verifier shares the #802 relation and this builder, this is deterministic replay and independent byte/basis accounting, not a second proof of the parent group-law circuit. The n3 toy test additionally evaluates copy/inverse/double/generic witnesses and compares both byte counts to a physically written #804 CNF.

Decision labels are `BOTH_WITHIN_FROZEN_CAPS` only if both complete DAGs and target-O byte counts fit, and `CAPACITY_CENSORED` if either arm reaches a frozen representation cap with a reproducible receipt. `STOP` means no admissible decision. Even `BOTH_WITHIN_FROZEN_CAPS` authorizes only later separately frozen target, solver, signed-column, and PDP-yield experiments. It is not such an experiment itself.

## Release and commands

The present PR runs only Python syntax, frozen hashes, independent basis replay, and n3 toy comparison in CI. `run.py` fails closed while `release_main_head` is null. After the prerequisite merges and independent review, make a new frozen commit pinning the exact merged `origin/main`, all source/input hashes, and unchanged or explicitly reviewed caps. The measured runner requires authenticated `gh pr view` read access; verify each merge commit is an ancestor of that exact main head. Do not rely on an unauthenticated hosted CI environment for the measured gate.

After that separate release commit and review, the authorized commands are:

    python3 research/notes/ecc2k130/m10_export_capacity_20260925/ci_replay.py
    python3 research/notes/ecc2k130/m10_export_capacity_20260925/run.py --out /private/tmp/ecc2k130-m10-capacity-producer
    python3 research/notes/ecc2k130/m10_export_capacity_20260925/verify.py --evidence /private/tmp/ecc2k130-m10-capacity-producer --out /private/tmp/ecc2k130-m10-capacity-independent.json

Archive the exact frozen source head, complete producer directory including STOP receipts, and independent replay in a separate outcome PR. Do not change the map, factor dimensions, circuit, caps, or verdict rule after seeing an outcome without a new review and freeze.
