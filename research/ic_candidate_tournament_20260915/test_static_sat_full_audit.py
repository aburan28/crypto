"""Exercise independent orbit-row replay on actual source-bound F5 evidence."""
import unittest

from audit_static_sat_full import (audited_status, hash_to_curve,
                                   orbit_columns, relation_row)
from identity import curve_record
from pathlib import Path
from oracle import Curve, InvalidEvidence, rank
from run_generic_exact_yield_audit import PANEL, load_evidence, record
from tournament import read


class FullSatAuditMathTests(unittest.TestCase):
    def test_next_paired_public_point_is_frozen_without_a_scalar(self):
        path=Path(__file__).parent/'goal_20260924/paired-fresh-n17a1/target-panel.json'
        panel=read(path)
        files=load_evidence(read(PANEL))
        report=record(files,'jobs/n17a1/f5/stdout.json')
        curve=Curve(report['fixture'])
        self.assertEqual(panel['curve_id'],curve_record(report['fixture'])
                         ['curve']['curve_id'])
        self.assertEqual(hash_to_curve(curve,b'ic-paired-target-v1',
                                       panel['target_input']['seed']),
                         (1,(52411,72106)))
        self.assertFalse(panel['target_input']['known_scalar_supplied'])
        self.assertIsNone(curve.mul((52411,72106),curve.r))

    def test_conflict_budget_is_censored_not_solver_failure(self):
        cms=dict(returncode=15,timed_out=False,
                 command=['cryptominisat5','--maxconfl','1000000'])
        self.assertEqual(audited_status('SOLVER_ERROR',cms,
                                        b's INDETERMINATE\n'),
                         'CONFLICT_BUDGET_INCONCLUSIVE')
        self.assertEqual(audited_status('SOLVER_ERROR',cms,
                                        b's CRASH\n'),'SOLVER_ERROR')

    def test_real_relations_reconstruct_full_rank_and_column_logs(self):
        files=load_evidence(read(PANEL))
        report=record(files,'jobs/n17a1/f5/stdout.json')
        curve=Curve(report['fixture'])
        base=[curve.decode(point) for point in report['factor_base']]
        columns,mapping=orbit_columns(curve,base)
        self.assertEqual(len(columns),29)
        rows=[]
        for relation in report['relations']:
            row,rhs=relation_row(curve,base,mapping,columns,
                                 relation['a'],relation['points'])
            rows.append(row)
            logs={tuple(map(int,item['point'])):int(item['log'])
                  for item in report['column_logs']}
            self.assertEqual(sum(value*logs[point]
                                 for value,point in zip(row,columns))%curve.r,rhs)
        self.assertEqual(rank(rows,len(columns),curve.r),29)

    def test_bad_scalar_does_not_make_a_row(self):
        files=load_evidence(read(PANEL))
        report=record(files,'jobs/n17a1/f5/stdout.json')
        curve=Curve(report['fixture'])
        base=[curve.decode(point) for point in report['factor_base']]
        columns,mapping=orbit_columns(curve,base)
        relation=report['relations'][0]
        with self.assertRaises(InvalidEvidence):
            relation_row(curve,base,mapping,columns,
                         (relation['a']+1)%curve.r,relation['points'])


if __name__=='__main__':
    unittest.main()
