import TreeSplit2.Data

/-!
# Round 2: the registered verdict on the recorded cells

`report` is the registered rule (`Rule.lean`, sha256 `ec513363…b1e3`);
`cells` is generated from `runs/registered/cells.jsonl`. Kernel `decide`.
-/

namespace TreeSplit2

theorem cells_count : cells.length = 102 := by decide

/-- Validity gate passes; every registered group is CLOSED. -/
theorem registered_report :
    report cells =
      { gate := true, t2Frob2 := .closed, t2Nbw := .closed, t2PolyLarge := .closed,
        zPoly := .closed, zNbw := .closed, zFrob2 := .closed } := by
  decide

end TreeSplit2
