import importlib.util
import itertools
import json
from pathlib import Path
import unittest

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('fiber_census',HERE/'independent_set.py')
census=importlib.util.module_from_spec(spec);spec.loader.exec_module(census)

class IndependentSetTests(unittest.TestCase):
    def test_every_graph_through_five_variables_matches_full_enumeration(self):
        for n in range(6):
            edges=list(itertools.combinations(range(n),2))
            for code in range(1<<len(edges)):
                graph=[0]*n
                for j,(a,b) in enumerate(edges):
                    if code&(1<<j):graph[a]|=1<<b;graph[b]|=1<<a
                expected=0
                for mask in range(1<<n):
                    if all(not(mask&(1<<a) and mask&(1<<b)) for a,b in edges if graph[a]&(1<<b)):
                        expected=census.choose(expected,mask)
                self.assertEqual(census.maximum_independent(graph)[0],expected)
    def test_polynomial_parity_can_remove_an_interaction(self):
        self.assertEqual(census.interaction_graph([[3,3,1,0]],2),[0,0])
        self.assertEqual(census.interaction_graph([[3],[3]],2),[2,1])
    def test_frozen_census_replays_and_every_restriction_is_affine(self):
        root=HERE/'linear_fiber_census_01'
        for name,digest in json.loads((root/'manifest.json').read_text())['files'].items():self.assertEqual(census.sha(root/name),digest)
        result=json.loads((root/'results.json').read_text())
        self.assertEqual(result['input_manifest_sha256'],census.sha(HERE/'resource_probe_02/manifest.json'))
        self.assertEqual(len(result['cells']),12)
        for row in result['cells']:
            graph=census.interaction_graph(row['polys'],row['n'])
            self.assertEqual(graph,row['graph'])
            mask,work=census.maximum_independent(graph)
            self.assertEqual(mask,row['independent_mask']);self.assertEqual(work,row['exact_search_work'])
            for p in row['polys']:
                self.assertTrue(all((m&mask).bit_count()<=1 for m in p))
            self.assertIsNone(row['solver_cost']);self.assertIsNone(row['speed_ratio'])

if __name__=='__main__':unittest.main()
