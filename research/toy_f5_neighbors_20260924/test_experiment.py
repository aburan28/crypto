import json
import gzip
import random
import unittest
from itertools import product
from pathlib import Path
from curves import Curve, O, certified_neighbors, semaev3, decomposition_value
from matrix import Ring, Echelon, measure, step
from experiment import encode, ground_truth, summarize, digest


class MatrixTests(unittest.TestCase):
    def test_limits(self):
        for n in [0, 7, 32]:
            with self.assertRaises(ValueError):
                Ring(n)

    def test_mobius_roundtrip(self):
        ring = Ring(6)
        randomizer = random.Random(2026)
        for _ in range(20):
            values = [randomizer.randrange(2) for _ in range(64)]
            f = ring.anf(values)
            self.assertEqual([ring.evaluate(f,a) for a in range(64)], values)

    def test_boolean_cancellation(self):
        ring = Ring(2)
        self.assertEqual(ring.product(ring.pack([3,1,2]), 1), ring.pack([1]))

    def test_frobenius_pruning(self):
        ring = Ring(3)
        f = ring.pack([3,4,0])
        a, _ = step(ring,[f],4,False)
        b, report = step(ring,[f],4,True)
        self.assertEqual(a,b)
        self.assertEqual(report['pruned_rows'],1)

    def test_unsound_naive_counterexample(self):
        ring = Ring(2)
        f = ring.pack([3,1,2])
        report = measure(ring,[f])
        self.assertEqual(report['roots'],[0])

    def test_zero_and_unit_ideals(self):
        ring = Ring(3)
        self.assertEqual(measure(ring,[])['roots'],list(range(8)))
        self.assertEqual(measure(ring,[1])['roots'],[])

    def test_random_rowspaces_at_every_degree(self):
        randomizer = random.Random(901)
        for n in range(2,7):
            ring = Ring(n)
            for _ in range(10):
                polynomials = [ring.pack(m for m in ring.monos if m.bit_count()<=3 and randomizer.randrange(3)==0) for _ in range(4)]
                for degree in range(2*n+1):
                    self.assertEqual(step(ring,polynomials,degree,False)[0],step(ring,polynomials,degree,True)[0])
                self.assertEqual(measure(ring,polynomials)['status'],'VERIFIED')


class CurveTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source, cls.h, cls.neighbors = certified_neighbors()

    def test_certificates(self):
        self.assertEqual([phi.target.b for phi,_ in self.neighbors],[93,225,92,224])
        for _, cert in self.neighbors:
            self.assertEqual(cert['dual_forward_checked_points'],288)
            self.assertEqual(cert['dual_reverse_checked_points'],288)
            self.assertEqual(cert['direction'],'descending')

    def test_group_law(self):
        c = self.source
        for p,q,r in product(self.h[:8],repeat=3):
            self.assertEqual(c.add(c.add(p,q),r),c.add(p,c.add(q,r)))
        for p in self.h:
            self.assertEqual(c.add(p,c.neg(p)),O)

    def test_semaev_and_lifting_exhaustive(self):
        reps = sorted({min(p,self.source.neg(p)) for p in self.h if p!=O and p[0]!=0})[:4]
        for c, mapping in [(self.source,lambda p:p)]+[(phi.target,phi) for phi,_ in self.neighbors]:
            support = [(mapping(p),mapping(self.source.neg(p))) for p in reps]
            for m in [2,3]:
                for target in self.h:
                    r, generators, slots = encode(c,support,mapping(target),m,False)
                    truth,_ = ground_truth(c,support,mapping(target),slots)
                    roots = {a for a in range(1<<r.n) if all(r.evaluate(f,a)==0 for f in generators)}
                    self.assertFalse(set(truth)-roots)
                    cr, cg, cs = encode(c,support,mapping(target),m,True)
                    canonical = {a for a in roots if slots[a]==tuple(sorted(slots[a]))}
                    self.assertEqual(canonical,{a for a in range(1<<cr.n) if all(cr.evaluate(f,a)==0 for f in cg)})

    def test_frozen_evidence(self):
        here = Path(__file__).resolve().parent
        raw = here/'results/raw.json.gz'
        if not raw.exists():
            self.skipTest('production evidence not generated yet')
        data = json.loads(gzip.decompress(raw.read_bytes()))
        self.assertEqual(summarize(data),json.loads((here/'results/summary.json').read_text()))
        self.assertEqual(data['contract_sha256'],digest(json.loads((here/'contract.json').read_text())))
        import hashlib
        for name, sha in data['source_sha256'].items():
            self.assertEqual(hashlib.sha256((here/name).read_bytes()).hexdigest(),sha)
        self.assertEqual(len(data['cases']),1280)
        for c in data['cases']:
            self.assertEqual(len(set(c['repetition_sha256'])),1)
            self.assertEqual([t['rowspace_sha256'] for t in c['traces']['f4']], [t['rowspace_sha256'] for t in c['traces']['f5']])


if __name__ == '__main__':
    unittest.main()
