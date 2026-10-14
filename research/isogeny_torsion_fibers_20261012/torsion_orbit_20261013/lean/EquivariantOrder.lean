import Init

/-!
The linear action on the two rational functions u=(y+1)/x and v=(y-1)/x
induced by translation by T, once the curve-level addition identity is
established, is (u,v) ↦ (-v,u-v). This theorem checks the order-three
consequence over the integers. The curve-level identity is not encoded here.
-/

def cubicAction (uv : Int × Int) : Int × Int :=
  (-uv.2, uv.1 - uv.2)

theorem cubicAction_order_three (u v : Int) :
    cubicAction (cubicAction (cubicAction (u, v))) = (u, v) := by
  simp [cubicAction]
  omega

#print axioms cubicAction_order_three
