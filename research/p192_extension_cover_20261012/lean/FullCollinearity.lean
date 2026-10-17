import Mathlib.Data.Finset.Card

/-!
If a complete census says that every four-point quartic line with at least
three selected points is already in the source family, any other such line
has at least two points outside the selected set. The complete finite-field
census is supplied as an explicit hypothesis, not encoded as an axiom.
-/

namespace P192ExtensionCover

theorem nonsource_quartic_line_has_two_unselected
    {P L : Type*} [DecidableEq P] [DecidableEq L]
    (selected : Finset P) (source : Finset L) (cut : L → Finset P)
    (hcomplete : ∀ line, 3 ≤ ((cut line) ∩ selected).card → line ∈ source)
    (line : L) (hfour : (cut line).card = 4)
    (hnot : line ∉ source) :
    2 ≤ ((cut line) \ selected).card := by
  have h_at_most_two : ((cut line) ∩ selected).card ≤ 2 := by
    by_contra h
    have h_three : 3 ≤ ((cut line) ∩ selected).card := by omega
    exact hnot (hcomplete line h_three)
  have hpartition := Finset.card_inter_add_card_sdiff (cut line) selected
  rw [hfour] at hpartition
  omega

end P192ExtensionCover

#print axioms P192ExtensionCover.nonsource_quartic_line_has_two_unselected
