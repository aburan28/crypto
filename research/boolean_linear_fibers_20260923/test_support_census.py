import importlib.util
import json
from pathlib import Path
import unittest

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('support_census',HERE/'support_census.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

class SupportTests(unittest.TestCase):
    def test_bounded_census_matches_direct_subsets_on_small_polynomials(self):
        for n in range(6):
            for seed in range(32):
                polys=[]
                for e in range(n+2):
                    polys.append([m for m in range(1<<n) if m.bit_count()<=2 and ((m+1)*(e+3)*(seed+5))%11<3])
                got=module.census(polys,n);self.assertTrue(got['complete'])
                best={}
                for mask in range(1<<n):
                    if any((mask&m).bit_count()==2 for p in polys for m in p):continue
                    untouched=sum(1<<e for e,p in enumerate(polys) if all(m&mask==0 for m in p))
                    key=mask.bit_count();value={'selected_mask':mask,'untouched_mask':untouched,'untouched_equations':untouched.bit_count()}
                    if key not in best or (value['untouched_equations'],-mask)>(best[key]['untouched_equations'],-best[key]['selected_mask']):best[key]=value
                self.assertEqual(got['best_by_size'],best)
    def test_caps_and_parity_do_not_create_an_optimum_claim(self):
        self.assertFalse(module.census([[]],6,1)['complete'])
        self.assertTrue(module.census([[]],0,1)['complete'])
        got=module.census([[3,3,1,1]],2)
        self.assertEqual(got['dependencies'],[0,0])
        self.assertEqual(got['best_by_size'][2]['untouched_equations'],1)
    def test_frozen_census_replays_without_a_performance_claim(self):
        root=HERE/'support_census_01'
        for name,digest in json.loads((root/'manifest.json').read_text())['files'].items():self.assertEqual(module.sha(root/name),digest)
        result=json.loads((root/'results.json').read_text())
        self.assertEqual(result['input_manifest_sha256'],module.sha(HERE/'resource_probe_02/manifest.json'))
        self.assertEqual(len(result['cells']),12)
        for c in result['cells']:
            got=module.census(c['polys'],c['n'],c['cap'])
            self.assertEqual(got['complete'],c['complete']);self.assertEqual(got['visited_sets'],c['visited_sets'])
            self.assertEqual({str(k):v for k,v in got['best_by_size'].items()},c['best_by_size'])
            self.assertIsNone(c['solver_cost']);self.assertIsNone(c['speed_ratio'])

if __name__=='__main__':unittest.main()
