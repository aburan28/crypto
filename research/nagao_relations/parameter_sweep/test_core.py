"""Post-run correctness audit; never imported by a timed solver."""
from collections import Counter
import itertools
import json
from pathlib import Path
import time
import unittest

import core

HERE=Path(__file__).resolve().parent
CONTRACT=json.loads((HERE/'contract.json').read_text())


class ExactSweepAudit(unittest.TestCase):
    def test_dynamic_program_counts_match_direct_signed_enumeration(self):
        p=CONTRACT['profiles'][0];truth=core.Truth(p);curve=truth.curve
        for ell,k in itertools.product(p['dimensions'],CONTRACT['summands']):
            base=core.factors(curve,core.support(truth.f,p,ell));actual=Counter()
            for chosen in itertools.combinations(base,k):
                xs={x for x,y in chosen}
                for signed in itertools.product(*[(point,curve.neg(point)) for point in chosen]):
                    total=None
                    for point in signed:total=curve.add(total,point)
                    if total is not None and total[0] not in xs:actual[total]+=1
            expected=truth.census(ell,k)[1]
            self.assertEqual({point:actual[point] for point in truth.points},expected)

    def test_subfield_projection_obstruction(self):
        for p in CONTRACT['profiles'][2:]:
            truth=core.Truth(p);f=truth.f;curve=truth.curve
            base=core.factors(curve,core.support(f,p,1))
            subfieldOrder=2*len(base)+2
            self.assertEqual(truth.h%subfieldOrder,0)
            for point in base:
                self.assertEqual(f.frob(point[1],3),point[1])
                self.assertIsNone(curve.mul(point,truth.h))
            self.assertGreater(truth.census(2,3)[0]['projected_columns'],0)

    def test_cold_solver_timeout_is_retained(self):
        p=CONTRACT['profiles'][0]
        target=json.loads((HERE/'targets.json').read_text())['profiles'][0]['targets'][0]['target']
        for variant in CONTRACT['stage_variants']:
            row=core.stageCell(p,3,3,target,variant,0)
            self.assertEqual(row['status'],'timeout')
            self.assertIsNone(row['witness'])
            self.assertGreater(row['field_api_counts']['totals']['multiplications'],0)

    def test_saved_full_dlp_scalars_on_actual_curves(self):
        truths={p['name']:core.Truth(p) for p in CONTRACT['profiles']};checks=0
        for line in (HERE/'results_v1/raw.jsonl').read_text().splitlines():
            row=json.loads(line)
            if row['kind'] not in ('dlp','rho') or row['status']!='verified':continue
            truth=truths[row['profile']];f=truth.f
            G=tuple(f.fromCoords(v) for v in row['generator']);Q=tuple(f.fromCoords(v) for v in row['target'])
            self.assertTrue(truth.curve.onCurve(G));self.assertTrue(truth.curve.onCurve(Q))
            self.assertIsNone(truth.curve.mul(G,row['subgroup_order']))
            self.assertEqual(truth.curve.mul(G,row['recovered_scalar']),Q);checks+=1
        self.assertEqual(checks,94)


if __name__=='__main__':unittest.main(verbosity=2)
