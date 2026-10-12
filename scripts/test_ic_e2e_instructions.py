"""Tests for scripts/ic_e2e_instructions.py (parsing and attribution; no valgrind needed)."""
from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import ic_e2e_instructions as ir  # noqa: E402

ANNOTATE = """--------------------------------------------------------------------------------
Profile data file 'x' (creator: callgrind-3.22.0)
Events recorded:  Ir
--------------------------------------------------------------------------------
Ir
6,654,044,047 (100.0%)  PROGRAM TOTALS

--------------------------------------------------------------------------------
Ir                      file:function
--------------------------------------------------------------------------------
3,272,031,189 (49.17%)  ???:<crypto_lib::cryptanalysis::koblitz_index_calculus::RelationCollector>::with_pair [ic]
  226,985,577 ( 3.41%)  ???:crypto_lib::cryptanalysis::koblitz_index_calculus::koblitz_signed_frobenius_rho_with_progress [ic]
  150,770,530 ( 2.27%)  ???:<alloc::vec::Vec<(FastPoint, u64, u64)> as SpecFromIterNested<...signed_rho_fast::{closure#1}>>::from_iter [ic]
    1,000,000 ( 0.02%)  ???:<crypto_lib::cryptanalysis::koblitz_index_calculus::PairSumTable>::build_within [ic]
"""


class ParseTests(unittest.TestCase):
    def test_total_from_summary_line_uses_the_ir_column(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "c.out"
            p.write_text("# callgrind format\nevents: Ir Dr\nfn=main\n0 5 1\nsummary: 123456 99\n")
            self.assertEqual(ir.parse_callgrind_total(p), 123456)
            p.write_text("events: Ir\ntotals: 42\n")
            self.assertEqual(ir.parse_callgrind_total(p), 42)
            p.write_text("events: Ir\n")
            with self.assertRaises(ValueError):
                ir.parse_callgrind_total(p)

    def test_inclusive_rows_and_symbol_lookup(self):
        rows = []
        for line in ANNOTATE.splitlines():
            m = ir.re.match(r"\s*([\d,]+)\s+\([\d. ]+%\)\s+\S+?:(.*)$", line)
            if m and "PROGRAM TOTALS" not in line:
                rows.append((int(m.group(1).replace(",", "")), m.group(2).strip()))
        self.assertEqual(len(rows), 4)
        self.assertEqual(ir.inclusive_of(rows, ir.RHO_SYMBOL), 226985577)
        self.assertEqual(ir.inclusive_of(rows, "RelationCollector>::with_pair"), 3272031189)
        self.assertEqual(ir.inclusive_of(rows, "PairSumTable>::build_within"), 1000000)
        self.assertIsNone(ir.inclusive_of(rows, "workflow::solve"))
        # the rho closure's from_iter row must not be mistaken for the rho routine
        self.assertNotEqual(ir.inclusive_of(rows, ir.RHO_SYMBOL), 150770530)

    def test_markdown_renders_one_row_per_rung(self):
        result = {"unit": "u", "rungs": [{"name": "a", "subgroup_bits": 20.5, "targets": 32, "ir_ic": 6414827561, "ir_rho": 226990004,
                                           "S_ic": 167088.0, "S_rho": 5912.4, "rho_over_ic_ir": 0.0354, "rho_ir_per_group_addition": 41604.0}]}
        text = ir.render_markdown(result)
        self.assertEqual(text.count("| `a` |"), 1)
        self.assertIn("6,414,827,561", text)
        self.assertIn("**0.035**", text)


if __name__ == "__main__":
    unittest.main()
