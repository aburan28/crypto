"""Check the proposed SAT relation matrix against real source-bound F5 rows."""
import unittest

from oracle import Curve, InvalidEvidence
from run_generic_exact_yield_audit import PANEL, load_evidence, record
from static_sat_matrix import RelationMatrix
from tournament import read


class StaticSatMatrixTests(unittest.TestCase):
    def test_real_f5_relations_recover_every_folded_column(self):
        files = load_evidence(read(PANEL))
        report = record(files, 'jobs/n17a1/f5/stdout.json')
        curve = Curve(report['fixture'])
        base = [curve.decode(point) for point in report['factor_base']]
        matrix = RelationMatrix(curve, base)
        self.assertEqual(len(base), 63)
        self.assertEqual(len(matrix.columns), 29)
        relations = report['relations']
        self.assertEqual(len(relations), 29)
        for relation in relations:
            matrix.push(relation['a'], relation['points'])
        self.assertEqual(matrix.rank, 29)
        logs = matrix.solve()
        reported = {tuple(map(int, entry['point'])): int(entry['log'])
                    for entry in report['column_logs']}
        self.assertEqual({point: log for point,log
                          in zip(matrix.columns, logs)}, reported)
        target_scalar = 2718
        target = curve.mul(curve.g, target_scalar)
        self.assertEqual(matrix.recover_target(
            logs, target, (63975-target_scalar) % curve.r, 1,
            (44, 58, 9)), target_scalar)

    def test_false_group_row_is_refused(self):
        files = load_evidence(read(PANEL))
        report = record(files, 'jobs/n17a1/f5/stdout.json')
        curve = Curve(report['fixture'])
        base = [curve.decode(point) for point in report['factor_base']]
        matrix = RelationMatrix(curve, base)
        first = report['relations'][0]
        with self.assertRaises(InvalidEvidence):
            matrix.push((first['a']+1) % curve.r, first['points'])
        self.assertEqual(matrix.rank, 0)


if __name__ == '__main__':
    unittest.main()
