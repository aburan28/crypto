-------------------------- MODULE Certification --------------------------
EXTENDS Naturals, FiniteSets, TLC
CONSTANT MaxOrder
VARIABLES lower, upper, actual, claimed, point
vars == <<lower,upper,actual,claimed,point>>
Init == /\ lower \in 1..MaxOrder
        /\ upper \in lower..MaxOrder
        /\ actual \in lower..upper
        /\ claimed \in lower..upper
        /\ point \in 0..(actual-1)
PointOrder == CHOOSE d \in 1..actual :
                 /\ (d*point) % actual = 0
                 /\ \A e \in 1..(d-1) : (e*point) % actual # 0
Annihilates == (claimed*point) % actual = 0
Certified == Annihilates /\ PointOrder > upper-lower
Sound == Certified => claimed = actual
SampleSound == Annihilates => claimed = actual
NontrivialSampleSound == (actual >= 2 /\ point > 0) => SampleSound
Next == UNCHANGED vars
Spec == Init /\ [][Next]_vars
==========================================================================
