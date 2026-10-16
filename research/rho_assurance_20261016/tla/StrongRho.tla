----------------------------- MODULE StrongRho -----------------------------
(***************************************************************************)
(* The strong distinguished-point Pollard rho of                           *)
(* src/cryptanalysis/koblitz_strong_rho.rs (StrongRhoG::solve), with the    *)
(* prime-order subgroup <G> replaced by the additive group Z_R:            *)
(*   G = 1, Q = d, a point is a residue x in 0..R-1, 0 is infinity,         *)
(*   k.P = k*P mod R, P + P' = P + P' mod R.                                *)
(* The signed-Frobenius automorphisms are x |-> +-LAM^k * x with LAM of     *)
(* multiplicative order NF mod R.  See README.md for the line-by-line map.  *)
(***************************************************************************)
EXTENDS Integers, Sequences, FiniteSets, TLC

CONSTANTS
    R,          \* prime subgroup order
    LAM,        \* Frobenius eigenvalue, multiplicative order NF mod R
    NF,         \* order of LAM (the field degree n in the Rust)
    L,          \* number of lockstep lanes
    NJ,         \* number of jumps J (2 or 3)
    J1, J2, J3, \* jump scalars s_j (jump point s_j.G); J3 unused when NJ = 2
    C1, C2,     \* partition(x) = (x*C1 + C2) % J
    C3, C4, D,  \* distinguished iff (x*C3 + C4) % D = 0
    WalkCap,    \* walk_cap (8 << dp_bits in the Rust)
    StepCap,    \* step_cap (global)
    ESCAPE,     \* params.escape_fruitless
    LaneOrderND,\* TRUE: lanes in a phase are taken in any order; FALSE: in index order (as the Rust)
    StrideLo, StrideHi, CursorLo, CursorHi, DLo, DHi, \* ranges explored by Init
    MUTANT      \* "none" or a deliberate bug, see below

ASSUME R > 2 /\ NF >= 1 /\ L >= 1 /\ D >= 1 /\ StepCap >= 1 /\ WalkCap >= 0
ASSUME NJ \in {2, 3}
JumpA == SubSeq(<<J1, J2, J3>>, 1, NJ)
ASSUME \A j \in 1..Len(JumpA) : JumpA[j] \in 1..R-1
ASSUME 1 <= StrideLo /\ StrideHi <= R-1 /\ 0 <= CursorLo /\ CursorHi <= R-1
ASSUME 0 <= DLo /\ DHi <= R-1
ASSUME MUTANT \in {"none", "canon_b", "sign", "sign_noverify", "escape_b"}

J      == Len(JumpA)
Pts    == 0..R-1          \* points (0 = infinity)
Scal   == 0..R-1          \* scalars mod R
Lanes  == 1..L
Md(x)  == ((x % R) + R) % R

\* lam_pow table and the automorphism multipliers {+-LAM^k}.
RECURSIVE Pow(_, _)
Pow(b, e) == IF e = 0 THEN 1 ELSE Md(Pow(b, e - 1) * b)
LamPow == [k \in 0..NF-1 |-> Pow(LAM, k)]
ASSUME Md(LamPow[NF-1] * LAM) = 1     \* "lambda must have order dividing n"
Mults  == {LamPow[k] : k \in 0..NF-1} \cup {Md(-LamPow[k]) : k \in 0..NF-1}
A      == Cardinality(Mults)          \* automorphisms

Min(S) == CHOOSE x \in S : \A y \in S : x <= y

\* canonicalize: least orbit member and the multiplier m with m*x = least.
\* Infinity (x = 0) is returned unchanged (multiplier 1), as in the Rust.
CanonMul == [x \in Pts |->
               IF x = 0 THEN 1
               ELSE LET best == Min({Md(m * x) : m \in Mults})
                    IN  CHOOSE m \in Mults : Md(m * x) = best]
Canon(x) == Md(CanonMul[x] * x)

CanonState(x, a, b) ==
    LET m == CanonMul[x] IN
    [pt |-> Md(m * x),
     a  |-> Md(a * m),
     b  |-> IF MUTANT = "canon_b" THEN b ELSE Md(b * m)]   \* MUTANT canon_b: b not multiplied

Inv == [v \in 1..R-1 |-> CHOOSE w \in 1..R-1 : Md(v * w) = 1]

Partition(x) == ((x * C1 + C2) % J) + 1          \* 1-based index into JumpA
IsDP(x)      == (x * C3 + C4) % D = 0

NoPt == 0   \* RawPointG::Infinity, also the sentinel that fills `previous`

LaneRec == [pt : Pts, a : Scal, b : Scal, prev : [1..4 -> Pts], nprev : 0..4, len : Nat, live : BOOLEAN]
IdleLane == [pt |-> NoPt, a |-> 0, b |-> 0, prev |-> <<NoPt, NoPt, NoPt, NoPt>>, nprev |-> 0, len |-> 0, live |-> FALSE]

VARIABLES
    d, strideA, cursorA, cursor,   \* target and walk-start generator
    lane, table,
    phase, pending,
    steps, walks, fruitless, escapes, capped, wasted, failed,
    sentinel,                      \* ghost: fruitless detections caused only by the Infinity sentinel
    longest,                       \* ghost: longest genuine cycle a fruitless detection has caught
    status, result

vars == <<d, strideA, cursorA, cursor, lane, table, phase, pending,
          steps, walks, fruitless, escapes, capped, wasted, failed, sentinel, longest, status, result>>
counters == <<steps, walks, fruitless, escapes, capped, wasted, failed, sentinel, longest>>

Init ==
    /\ d \in DLo..DHi
    /\ strideA \in StrideLo..StrideHi            \* stride_a in [1, r)
    /\ cursorA \in CursorLo..CursorHi            \* cursor_a in [0, r)
    /\ cursor = Md(cursorA * 1 + d)              \* cursor = cursor_a.G + Q
    /\ lane = [l \in Lanes |-> IdleLane]
    /\ table = <<>>                              \* empty function
    /\ phase = "begin"
    /\ pending = {}
    /\ steps = 0 /\ walks = 0 /\ fruitless = 0 /\ escapes = 0
    /\ capped = 0 /\ wasted = 0 /\ failed = 0 /\ sentinel = 0 /\ longest = 0
    /\ status = "running"
    /\ result = 0

\* Which pending lane may go next.
Pick(l) == /\ l \in pending
           /\ LaneOrderND \/ l = Min(pending)

Running == status = "running"

(* while steps < step_cap *)
Begin ==
    /\ Running /\ phase = "begin"
    /\ IF steps >= StepCap
         THEN /\ status' = "capped"
              /\ UNCHANGED <<phase, pending>>
         ELSE /\ phase' = "start"
              /\ pending' = {l \in Lanes : ~lane[l].live}
              /\ UNCHANGED status
    /\ UNCHANGED <<d, strideA, cursorA, cursor, lane, table, counters, result>>

(* Phase 1: one turn of the start loop for an idle lane. *)
StartLane(l) ==
    /\ Running /\ phase = "start" /\ Pick(l)
    /\ walks' = walks + 1
    /\ cursor' = Md(cursor + strideA)            \* cursor = add(cursor, stride)
    /\ cursorA' = Md(cursorA + strideA)
    /\ IF cursor = NoPt
         THEN UNCHANGED <<lane, pending>>        \* `continue`: same lane, next start
         ELSE LET c == CanonState(cursor, cursorA, 1) IN
              /\ lane' = [lane EXCEPT ![l] = [pt |-> c.pt, a |-> c.a, b |-> c.b,
                                             prev |-> <<NoPt, NoPt, NoPt, NoPt>>,
                                             nprev |-> 0,
                                             len |-> 0, live |-> TRUE]]
              /\ pending' = pending \ {l}
    /\ UNCHANGED <<d, strideA, table, phase, steps, fruitless, escapes, capped,
                   wasted, failed, sentinel, longest, status, result>>

EndStart ==
    /\ Running /\ phase = "start" /\ pending = {}
    /\ phase' = "dp" /\ pending' = Lanes
    /\ UNCHANGED <<d, strideA, cursorA, cursor, lane, table, counters, status, result>>

(* Phase 2: a lane standing on a distinguished point. *)
DpLane(l) ==
    /\ Running /\ phase = "dp" /\ Pick(l)
    /\ pending' = pending \ {l}
    /\ LET s == lane[l] IN
       IF ~(s.live /\ IsDP(s.pt))
         THEN UNCHANGED <<lane, table, wasted, failed, status, result>>
         ELSE
         /\ lane' = [lane EXCEPT ![l].live = FALSE]
         /\ IF s.pt \notin DOMAIN table
              THEN /\ table' = table @@ (s.pt :> <<s.a, s.b>>)
                   /\ UNCHANGED <<wasted, failed, status, result>>
              ELSE
              LET h   == table[s.pt]
                  den == Md(s.b - h[2])
              IN
              /\ UNCHANGED table
              /\ IF den = 0
                   THEN /\ wasted' = wasted + 1
                        /\ UNCHANGED <<failed, status, result>>
                   ELSE
                   LET cand == IF MUTANT \in {"sign", "sign_noverify"}
                                 THEN Md((s.a - h[1]) * Inv[den])   \* MUTANT: wrong sign
                                 ELSE Md((h[1] - s.a) * Inv[den])
                   IN
                   IF Md(cand * 1) = d \/ MUTANT = "sign_noverify"  \* scalar_mul(G, cand) == q
                     THEN /\ status' = "solved" /\ result' = cand
                          /\ UNCHANGED <<wasted, failed>>
                     ELSE /\ failed' = failed + 1
                          /\ UNCHANGED <<wasted, status, result>>
    /\ UNCHANGED <<d, strideA, cursorA, cursor, phase, steps, walks, fruitless, escapes, capped, sentinel, longest>>

EndDp ==
    /\ Running /\ phase = "dp" /\ pending = {}
    /\ phase' = "step" /\ pending' = {l \in Lanes : lane[l].live}
    /\ UNCHANGED <<d, strideA, cursorA, cursor, lane, table, counters, status, result>>

(* Phase 3: one step of a live lane. *)
StepLane(l) ==
    /\ Running /\ phase = "step" /\ Pick(l)
    /\ pending' = pending \ {l}
    /\ steps' = steps + 1
    /\ LET s    == lane[l]
           s_j  == JumpA[Partition(s.pt)]
           nxt  == CanonState(Md(s.pt + s_j), Md(s.a + s_j), s.b)
           len1 == s.len + 1
           prevSet == {s.prev[i] : i \in 1..4}
           realPrev == {s.prev[i] : i \in 1..s.nprev}   \* ghost: slots really written
           cycle == nxt.pt = s.pt \/ nxt.pt \in prevSet
       IN
       /\ sentinel' = IF cycle /\ nxt.pt # s.pt /\ nxt.pt \notin realPrev
                        THEN sentinel + 1 ELSE sentinel
       /\ longest' = LET cyc == IF nxt.pt = s.pt THEN 1
                                ELSE IF \E i \in 1..s.nprev : s.prev[i] = nxt.pt
                                     THEN 1 + Min({i \in 1..s.nprev : s.prev[i] = nxt.pt})
                                     ELSE 0
                       IN IF cycle /\ cyc > longest THEN cyc ELSE longest
       /\ IF cycle
            THEN IF ~ESCAPE
                   THEN /\ fruitless' = fruitless + 1
                        /\ lane' = [lane EXCEPT ![l].len = len1, ![l].live = FALSE]
                        /\ UNCHANGED <<escapes, capped>>
                   ELSE LET e == CanonState(Md(2 * nxt.pt), Md(2 * nxt.a),
                                            IF MUTANT = "escape_b" THEN nxt.b ELSE Md(2 * nxt.b))
                        IN /\ escapes' = escapes + 1
                           /\ lane' = [lane EXCEPT ![l] = [pt |-> e.pt, a |-> e.a, b |-> e.b,
                                                          prev |-> <<NoPt, NoPt, NoPt, NoPt>>,
                                                          nprev |-> 0,
                                                          len |-> len1, live |-> TRUE]]
                           /\ UNCHANGED <<fruitless, capped>>
          ELSE IF len1 > WalkCap
            THEN /\ capped' = capped + 1
                 /\ lane' = [lane EXCEPT ![l].len = len1, ![l].live = FALSE]
                 /\ UNCHANGED <<fruitless, escapes>>
          ELSE /\ lane' = [lane EXCEPT ![l] = [pt |-> nxt.pt, a |-> nxt.a, b |-> nxt.b,
                                              prev |-> <<s.pt, s.prev[1], s.prev[2], s.prev[3]>>,
                                              nprev |-> IF s.nprev < 4 THEN s.nprev + 1 ELSE 4,
                                              len |-> len1, live |-> TRUE]]
               /\ UNCHANGED <<fruitless, escapes, capped>>
    /\ UNCHANGED <<d, strideA, cursorA, cursor, table, phase, walks, wasted, failed, status, result>>

EndStep ==
    /\ Running /\ phase = "step" /\ pending = {}
    /\ phase' = "begin" /\ pending' = {}
    /\ UNCHANGED <<d, strideA, cursorA, cursor, lane, table, counters, status, result>>

\* Explicit terminal state: solved (Some) or capped (None).
Done == ~Running /\ UNCHANGED vars

Next ==
    \/ Begin \/ EndStart \/ EndDp \/ EndStep
    \/ \E l \in Lanes : StartLane(l) \/ DpLane(l) \/ StepLane(l)
    \/ Done

Spec     == Init /\ [][Next]_vars
FairSpec == Spec /\ WF_vars(Next)

-----------------------------------------------------------------------------
(* Properties *)

TypeOK ==
    /\ d \in Scal /\ strideA \in 1..R-1 /\ cursorA \in Scal /\ cursor \in Pts
    /\ lane \in [Lanes -> LaneRec]
    /\ DOMAIN table \subseteq Pts
    /\ \A k \in DOMAIN table : table[k] \in Scal \X Scal
    /\ phase \in {"begin", "start", "dp", "step"}
    /\ pending \subseteq Lanes
    /\ {steps, walks, fruitless, escapes, capped, wasted, failed} \subseteq Nat
    /\ status \in {"running", "solved", "capped"}
    /\ result \in Scal

\* point = a.G + b.Q for every lane state (live or not) and every table entry,
\* and the walk-start cursor = cursor_a.G + Q.
RepInv ==
    /\ \A l \in Lanes : lane[l].pt = Md(lane[l].a + lane[l].b * d)
    /\ \A k \in DOMAIN table : k = Md(table[k][1] + table[k][2] * d)
    /\ cursor = Md(cursorA + d)

\* Every live lane stands on the canonical representative of its orbit;
\* table keys are canonical distinguished points.
CanonInv ==
    /\ \A l \in Lanes : lane[l].live => Canon(lane[l].pt) = lane[l].pt
    /\ \A k \in DOMAIN table : Canon(k) = k /\ IsDP(k)

NoFailedCollision == failed = 0

SolvedCorrect == status = "solved" => result = d

\* Diagnostic (not part of correctness): a live lane never exceeds walk_cap.
WalkCapRespected == \A l \in Lanes : lane[l].live => lane[l].len <= WalkCap

\* Liveness under weak fairness.
Terminates == <>(status # "running")

-----------------------------------------------------------------------------
(* Reachability witnesses.  Each is an invariant that TLC is EXPECTED to     *)
(* violate: the counterexample is a trace that exercises the named branch,  *)
(* so the safety properties above are not vacuously true.                   *)

W_Solved        == status # "solved"
W_StepCapped    == status # "capped"
W_Wasted        == wasted = 0                  \* table hit with b = b'
W_Fruitless     == fruitless = 0               \* walk abandoned in a fruitless cycle
W_Escape        == escapes = 0                 \* fruitless cycle escaped by doubling
W_WalkCapped    == capped = 0                  \* walk abandoned at walk_cap
W_InfinityStart == ~(phase = "start" /\ pending # {} /\ cursor = NoPt)  \* start loop `continue`
W_LaneAtInfinity == \A l \in Lanes : lane[l].live => lane[l].pt # NoPt
W_Sentinel      == sentinel = 0                \* fruitless only via the Infinity sentinel
W_TableTwo      == Cardinality(DOMAIN table) < 2
W_Cycle5        == longest < 5                 \* a genuine 5-cycle is caught as fruitless
=============================================================================
