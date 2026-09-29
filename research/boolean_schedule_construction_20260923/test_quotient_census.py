import importlib.util
import json
from pathlib import Path
import unittest

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('quotient_census',HERE/'quotient_census.py')
q=importlib.util.module_from_spec(spec);spec.loader.exec_module(q)
def value(coefficients,point):
    result=0
    for m,c in coefficients.items():
        if point&m==m:result^=c
    return result

class QuotientTests(unittest.TestCase):
    def test_every_small_projection_is_linear_and_annihilates_its_span(self):
        for code in range(256):
            columns=[v for v in range(8) if code&(1<<v)];rows=q.basis(columns)
            for c in columns:self.assertEqual(q.project(c,rows),0)
            for a in range(8):
                for b in range(8):self.assertEqual(q.project(a^b,rows),q.project(a,rows)^q.project(b,rows))
    def test_projection_preserves_necessity_but_not_sufficiency(self):
        # xy+1 projects to zero; only x=y=1 solves the original equation.
        c=q.analyze([[3,0]],2,2)
        self.assertEqual(c['projected_coefficients'],{})
        original=[x for x in range(4) if value(c['original_coefficients'],x)==0]
        projected=[x for x in range(4) if value(c['projected_coefficients'],x)==0]
        self.assertEqual(original,[3]);self.assertEqual(projected,[0,1,2,3])
    def test_all_small_polynomials_have_exact_projected_evaluations(self):
        monomials=[0,1,2,4,3,5,6]
        for code in range(128):
            first=[m for i,m in enumerate(monomials) if code&(1<<i)]
            for second in [[],[0],[1,0],[3,5,2],[6,1,4]]:
                for k in [1,2,3]:
                    c=q.analyze([first,second],3,k);rows=c['elimination_basis']
                    for x in range(8):
                        expected=q.project(value(c['original_coefficients'],x),rows)
                        self.assertEqual(value(c['projected_coefficients'],x),expected)
                        if value(c['original_coefficients'],x)==0:self.assertEqual(expected,0)
    def test_frozen_discovery_census_replays(self):
        root=HERE/'quotient_census_01'
        for name,digest in json.loads((root/'manifest.json').read_text())['files'].items():self.assertEqual(q.sha(root/name),digest)
        result=json.loads((root/'results.json').read_text());self.assertEqual(len(result['cells']),72)
        self.assertEqual(result['input_manifest_sha256'],q.sha(HERE/'run_01/manifest.json'))
        for cell in result['cells']:
            source=HERE/'run_01'/f"n{cell['n']}-{cell['seed']}-{cell['family']}.jsonl"
            self.assertEqual(q.sha(source),cell['input_sha256'])
            with source.open() as f:fixture=json.loads(next(f))
            expected=q.analyze(fixture['polys'],cell['n'],cell['low_variables'])
            for key,value in expected.items():self.assertEqual(json.loads(json.dumps(value)),cell[key])
            self.assertIsNone(cell['performance_ratio']);self.assertFalse(cell['complete_solver_implemented'])

if __name__=='__main__':unittest.main()
