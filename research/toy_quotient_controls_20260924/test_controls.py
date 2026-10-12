import json
import random
import unittest
from itertools import combinations_with_replacement, product
from pathlib import Path
import algebra as a
from curves import certified_neighbors, O, decomposition_value, MUL, SQUARE
from independent import audit
from run import pipeline, oracle, opcode_call

CONTRACT = json.loads((Path(__file__).resolve().parent/'contract.json').read_text())


class AlgebraTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source, cls.h, cls.neighbors = certified_neighbors()
        choices = sorted({min(p,cls.source.neg(p)) for p in cls.h if p!=O and p[0]!=0})
        cls.support = [(p,cls.source.neg(p)) for p in random.Random(101).sample(choices,4)]

    def test_symmetric_equation_identity(self):
        rng = random.Random(31415)
        for b in (1,93,225,92,224):
            for m in (2,3):
                for _ in range(150):
                    xs = [rng.randrange(256) for _ in range(m)]
                    target = O if _%7==0 else (rng.randrange(256),0)
                    self.assertEqual(a.symmetric_value(a.elementary(xs),target,b),decomposition_value(xs,target,b))
                for x in range(256):
                    xs = [x]*m
                    self.assertEqual(a.symmetric_value(a.elementary(xs),O,b),decomposition_value(xs,O,b))

    def test_repeated_root_recovery(self):
        for m in (2,3):
            seen = set()
            for indices in combinations_with_replacement(range(4),m):
                e = a.elementary([self.support[i][0][0] for i in indices])
                self.assertNotIn(e,seen)
                seen.add(e)
                self.assertEqual(a.recover(e,self.support),indices)
            self.assertEqual(len(seen),10 if m==2 else 20)

    def test_invalid_support_and_non_split_polynomial(self):
        with self.assertRaises(ValueError):
            a.setup(self.support*2,3,'invariant',CONTRACT)
        with self.assertRaises(ValueError):
            a.recover((0,0),self.support)  # T^2, support excludes zero

    def test_coordinate_scaling_identity(self):
        rng = random.Random(119)
        for b in (1,93,225,92,224):
            for m in (2,3):
                for _ in range(100):
                    xs = [rng.randrange(256) for _ in range(m)]
                    target = O if _%5==0 else (rng.randrange(256),0)
                    exponent = (2 if m==2 else 8) if target==O else (8 if m==2 else 24)
                    self.assertEqual(a.scaled_value(xs,target,b,2),MUL[a.power(2,exponent)][decomposition_value(xs,target,b)])

    def test_scaled_point_equations(self):
        u=2
        for curve,points in [(self.source,self.h)]+[(phi.target,[phi(p) for p in self.h]) for phi,_ in self.neighbors]:
            for p in points:
                if p==O:
                    continue
                x,y=MUL[a.power(u,2)][p[0]],MUL[a.power(u,3)][p[1]]
                self.assertEqual(SQUARE[y]^MUL[u][MUL[x][y]],MUL[SQUARE[x]][x]^MUL[a.power(u,6)][curve.b])

    def test_padding_excluded(self):
        for m in (2,3):
            setup=a.setup(self.support,m,'invariant',CONTRACT)
            ring,gens=a.encode(self.support,O,1,m,'invariant',setup,CONTRACT)
            self.assertEqual(ring.n,4 if m==2 else 5)
            for index in range(len(setup['domain']),1<<ring.n):
                self.assertEqual(ring.evaluate(gens[-1],index),1)

    def test_all_variants_against_point_oracle(self):
        for m in (2,3):
            for target in (O,self.h[7],self.h[31]):
                truth=oracle(self.support,target,m)
                degrees=[]
                for variant in CONTRACT['variants']:
                    result,_,_=pipeline(self.support,target,1,m,variant,CONTRACT)
                    self.assertEqual(result['solutions'],truth)
                    degrees.append(result['completion_degree'])
                self.assertEqual(degrees[0],degrees[4])  # generator permutation preserves V_D
                self.assertEqual(degrees[0],degrees[2])  # four-label permutation is affine over F2^2

    def test_opcode_replay_and_hook_cleanup(self):
        import sys
        first,count,_=pipeline(self.support,O,1,2,'invariant',CONTRACT,True)
        second,count2,_=pipeline(self.support,O,1,2,'invariant',CONTRACT,True)
        self.assertEqual(first,second)
        self.assertEqual(count,count2)
        self.assertTrue(all(value>0 for value in count.values()))
        def fails():
            raise ValueError('expected')
        with self.assertRaises(ValueError):
            opcode_call(fails)
        self.assertIsNone(sys.gettrace())

    def test_independent_f5b_and_buchberger(self):
        for variant in ('ordered','canonical','invariant'):
            for m in (2,3):
                result,_,_=pipeline(self.support,self.h[7],1,m,variant,CONTRACT)
                receipt=audit(result['n'],result['generators'],result['roots'])
                self.assertEqual(receipt['status'],'VERIFIED')
        for generators,roots in [([],list(range(4))),([1],[])]:
            self.assertEqual(audit(2,generators,roots)['status'],'VERIFIED')


    def test_independent_budget_and_negative_control(self):
        result,_,_=pipeline(self.support,self.h[7],1,2,'ordered',CONTRACT)
        receipt=audit(result['n'],result['generators'],result['roots'],call_cap=1)
        self.assertEqual(receipt['status'],'VERIFIED')
        self.assertEqual(receipt['f5b_status'],'BUDGET_EXCEEDED')
        with self.assertRaises(AssertionError):
            audit(2,[],[])


if __name__=='__main__':
    unittest.main()
