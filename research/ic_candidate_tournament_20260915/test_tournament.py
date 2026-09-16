import copy
import json
from pathlib import Path
import tempfile
import unittest

from oracle import Curve, InvalidEvidence, rank, verify
from tournament import comparison, gate, parse_profiles


class ArithmeticTests(unittest.TestCase):
    def setUp(self):
        # Independent tiny field and curve, with no dependency on Rust outputs.
        self.c = object.__new__(Curve)
        self.c.n = 5
        self.c.a = 0
        self.c.modulus = 0b100101

    def test_field_and_group_laws_exhaustively(self):
        c=self.c
        for x in range(1,32):
            self.assertEqual(c.fm(x,c.inv(x)),1)
        points=[None]
        for x in range(32):
            for y in range(32):
                try:points.append(c.decode([x,y]))
                except InvalidEvidence:pass
        for p in points:
            self.assertIsNone(c.add(p,c.neg(p)))
            for q in points:
                s=c.add(p,q)
                self.assertEqual(s,c.add(q,p))
                self.assertEqual(c.add(s,c.neg(q)),p)
                if s is not None:self.assertEqual(c.decode(list(s)),s)

    def test_rank_uses_subgroup_field_and_rejects_dependence(self):
        self.assertEqual(rank([[1,1],[2,2],[3,3]],2,7),1)
        self.assertEqual(rank([[1,1],[1,2]],2,7),2)
        self.assertEqual(rank([[1,1],[1,3]],2,2),1)


class ProfilingTests(unittest.TestCase):
    def write_profile(self, root, name, part, phase, total):
        (root/name).write_text(f'part: {part}\ndesc: Trigger: {phase}\nevents: Ir\nsummary: {total}\ntotals: {total}\n')

    def test_intervals_include_startup_and_cleanup_exactly_once(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            self.write_profile(root,'callgrind.out.1',1,'Client Request: startup_and_input',17)
            self.write_profile(root,'callgrind.out',2,'Program termination',11)
            (root/'stderr.txt').write_text('==9== Collected : 28\n')
            self.assertEqual(sum(parse_profiles(root).values()),28)
            (root/'stderr.txt').write_text('==9== Collected : 29\n')
            with self.assertRaises(InvalidEvidence):parse_profiles(root)

    def test_missing_final_interval_never_becomes_a_low_cost_win(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            self.write_profile(root,'callgrind.out.1',1,'Client Request: startup_and_input',17)
            (root/'stderr.txt').write_text('==9== Collected : 17\n')
            with self.assertRaises(InvalidEvidence):parse_profiles(root)


def rows(ratio=1):
    records=[]
    for c in range(4):
        for case in range(4):
            for rep in range(3):
                for arm,factor in [('incumbent',1),('candidate',ratio)]:
                    records.append({'arm':arm,'case':f'{c}-{case}','cell':str(c),
                        'repetition':rep,'case_sha256':f'{c}-{case}', 'status':'VERIFIED',
                        'total_operations':1000*factor,'native_process':{'process_wall_seconds':.1},
                        'certificate':{'factor_base_sha256':'fixed-base'}})
    return records


class DecisionTests(unittest.TestCase):
    def setUp(self):
        self.contract={'confirmation_ratio':.8,'max_cell_ratio':1.1}

    def test_aa_control_cannot_promote(self):
        result=comparison(rows(),'candidate',draws=200)
        self.assertEqual(result['candidate_over_baseline'],1)
        self.assertFalse(gate(result,self.contract))
        self.assertEqual(result['paired_cases'],16)  # Not 48 timing repetitions.
        self.assertEqual(result['independent_curve_blocks'],4)

    def test_qualified_cost_reduction(self):
        self.assertTrue(gate(comparison(rows(.7),'candidate',draws=200),self.contract))

    def test_timeout_and_missing_cost_block_promotion(self):
        for field,value in [('status','TIMEOUT'),('total_operations',None)]:
            data=rows(.5)
            data[1][field]=value
            self.assertFalse(gate(comparison(data,'candidate',draws=200),self.contract))

    def test_changed_target_rejected(self):
        data=rows(.5);data[1]['case_sha256']='different'
        with self.assertRaises(InvalidEvidence):comparison(data,'candidate',draws=200)

    def test_changed_factor_base_is_not_a_matched_implementation_gain(self):
        data=rows(.5)
        data[1]['certificate']['factor_base_sha256']='different-base'
        with self.assertRaises(InvalidEvidence):comparison(data,'candidate',draws=200)

    def test_per_cell_regression_blocks_aggregate_win(self):
        data=rows(.1)
        for row in data:
            if row['arm']=='candidate' and row['cell']=='0':row['total_operations']=1200
        result=comparison(data,'candidate',draws=200)
        self.assertLess(result['candidate_over_baseline'],.8)
        self.assertFalse(gate(result,self.contract))


if __name__=='__main__':unittest.main()
