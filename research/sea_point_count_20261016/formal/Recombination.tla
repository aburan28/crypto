------------------------- MODULE Recombination -------------------------
EXTENDS Naturals, Integers, FiniteSets, TLC

CONSTANTS MaxGroup, MaxCount, MaxRounds, Legacy
VARIABLES n, count, truth, point, stride, allowed, phase, j, baby, giant,
          hits, survivors, round, result
vars == <<n,count,truth,point,stride,allowed,phase,j,baby,giant,hits,
          survivors,round,result>>

Indices == 0..(count-1)
(* Exact integer ceil(sqrt(count)); all models have count > 0. *)
B == CHOOSE b \in 1..count : /\ b*b >= count /\ (b-1)*(b-1) < count
Step == (stride*point) % n
Base == (truth*stride*point) % n
Key(k) == (k*Step) % n
Expected == {k \in Indices : Key(k) = Base}
Kernel(p) == {k \in Indices : ((truth*stride + n - k*stride)*p) % n = 0}
(* Use nonnegative residues: subtraction may be negative; % is Euclidean. *)
GiantKey(i) == (Base - (i*B*Step)) % n

Init ==
  /\ n \in 1..MaxGroup
  /\ count \in 1..MaxCount
  /\ truth \in 0..(count-1)
  /\ point \in 0..(n-1)
  /\ stride \in {1,2,3}
  /\ allowed \in {0..(count-1), {truth},
                  {k \in 0..(count-1) : k % 2 = truth % 2}}
  /\ phase = "Baby"
  /\ j = 0
  /\ baby = [r \in 0..(n-1) |-> {}]
  /\ giant = 0
  /\ hits = {}
  /\ survivors = 0..(count-1)
  /\ round = 0
  /\ result = -1

BabyStep ==
  /\ phase = "Baby" /\ j < B
  /\ baby' = [baby EXCEPT ![Key(j)] =
        IF Legacy /\ @ # {} THEN @ ELSE @ \cup {j}]
  /\ j' = j+1
  /\ UNCHANGED <<n,count,truth,point,stride,allowed,phase,giant,hits,
                  survivors,round,result>>
BeginGiants ==
  /\ phase = "Baby" /\ j = B
  /\ phase' = "Giant"
  /\ UNCHANGED <<n,count,truth,point,stride,allowed,j,baby,giant,hits,
                  survivors,round,result>>
GiantStep ==
  /\ phase = "Giant" /\ giant*B < count
  /\ hits' = hits \cup {giant*B+b : b \in
                   {v \in baby[GiantKey(giant)] : giant*B+v < count}}
  /\ giant' = giant+1
  /\ UNCHANGED <<n,count,truth,point,stride,allowed,phase,j,baby,
                  survivors,round,result>>
FinishSearch ==
  /\ phase = "Giant" /\ giant*B >= count
  /\ survivors' = survivors \cap hits \cap allowed
  /\ phase' = "Choose"
  /\ UNCHANGED <<n,count,truth,point,stride,allowed,j,baby,giant,hits,round,result>>
ReturnUnique ==
  /\ phase = "Choose" /\ Cardinality(survivors) = 1
  /\ result' = CHOOSE k \in survivors : TRUE
  /\ phase' = "Done"
  /\ UNCHANGED <<n,count,truth,point,stride,allowed,j,baby,giant,hits,survivors,round>>
ReturnAmbiguous ==
  /\ phase = "Choose" /\ Cardinality(survivors) # 1 /\ round = MaxRounds
  /\ phase' = "Ambiguous"
  /\ UNCHANGED <<n,count,truth,point,stride,allowed,j,baby,giant,hits,survivors,round,result>>
NextSample ==
  /\ phase = "Choose" /\ Cardinality(survivors) # 1 /\ round < MaxRounds
  /\ \E p \in 0..(n-1) :
       \/ /\ survivors' = survivors \cap Kernel(p)
          /\ phase' = "Choose"
          /\ UNCHANGED <<point,j,baby,giant,hits>>
       \/ /\ point' = p
          /\ phase' = "Baby"
          /\ j' = 0 /\ giant' = 0 /\ hits' = {}
          /\ baby' = [r \in 0..(n-1) |-> {}]
          /\ UNCHANGED survivors
  /\ round' = round+1
  /\ UNCHANGED <<n,count,truth,stride,allowed,result>>
Terminal == phase \in {"Done","Ambiguous"}
Next == BabyStep \/ BeginGiants \/ GiantStep \/ FinishSearch \/
        ReturnUnique \/ ReturnAmbiguous \/ NextSample
Spec == Init /\ [][Next]_vars /\ WF_vars(Next)

TypeOK == /\ survivors \subseteq Indices /\ hits \subseteq Indices
          /\ phase \in {"Baby","Giant","Choose","Done","Ambiguous"}
TrueSurvives == truth \in survivors
NontrivialTrueSurvives == (n >= 2 /\ point > 0) => truth \in survivors
UniqueSound == phase = "Done" => result = truth
BabyComplete == phase = "Giant" =>
                  \A r \in 0..(n-1) : baby[r] = {v \in 0..(B-1) : Key(v)=r}
SearchPrefix == phase = "Giant" =>
                hits = {k \in Expected : k < giant*B}
Terminates == <>Terminal
=======================================================================
