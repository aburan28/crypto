"""Independent wide circuit witnesses; planted checks are not relation yield."""
import random
import unittest
import indexcalc_e2e as e


class WideWitnessTests(unittest.TestCase):
    def test_chained_s3_against_field_arithmetic(self):
        f, c, q, g, eigen, reps, lookup, prog, roots = e.setup(131, 2, 3, e.Ledger())
        rng = random.Random(91317)
        base = list(lookup)
        wide = False
        for _ in range(12):
            while True:
                points = rng.sample(base, 3)
                intermediate = c.add(points[0], points[1])
                target = c.add(intermediate, points[2])
                if intermediate is not None and target is not None and e.proper(points, [1]*3, c):
                    break
            coordinates = [f.toCoords(p[0]) for p in points]
            wide |= any(x >> 64 for x in coordinates)
            self.assertIsNotNone(e.strictLift(f, c, coordinates, target, lookup))
            values = dict(zip(('p0', 'p1', 'p2'), [p[0] for p in points]))
            values.update(t0=intermediate[0], r=target[0], one=f.one())
            for perturb in (False, True):
                if perturb:
                    values['r'] = f.fromCoords(rng.getrandbits(131))
                inputs = {(name, bit): (f.toCoords(value) >> bit) & 1
                          for name, value in values.items() for bit in range(131)}
                actual = prog.evaluate(inputs, roots)
                expected = []
                for u, v, w in ((values['p0'], values['p1'], values['t0']),
                                (values['t0'], values['p2'], values['r'])):
                    uv = f.mul(u, v)
                    s = f.add(f.add(uv, f.mul(u, w)), f.mul(v, w))
                    result = f.add(f.add(f.frob(s, 1), f.mul(uv, w)), f.one())
                    bits = f.toCoords(result)
                    expected.extend((bits >> bit) & 1 for bit in range(131))
                self.assertEqual(actual, expected)
                if not perturb:
                    self.assertFalse(any(actual))
                else:
                    self.assertTrue(any(actual))
        self.assertTrue(wide)


if __name__ == '__main__':
    unittest.main()
