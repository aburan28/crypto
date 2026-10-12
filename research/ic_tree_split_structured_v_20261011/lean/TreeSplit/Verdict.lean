import TreeSplit.Data

/-!
# The registered verdict on the recorded cells

`report` is the registered rule (`Rule.lean`, sha256 `714601c1…760d`);
`cells` is generated from `runs/registered/cells.jsonl`. The kernel
evaluates the rule on the records (`decide`), so `RESULTS.md` cannot state
a verdict the registered rule does not give on these data.
-/

namespace TreeSplit

theorem cells_count : cells.length = 64 := by decide

/-- Validity gate passes; `poly`, `frob` (`n = 21`) and `rand` are all CLOSED. -/
theorem registered_report :
    report cells = { gate := true, poly := .closed, frob := .closed, rand := .closed } := by
  decide

end TreeSplit
