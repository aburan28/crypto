"""Exhaustive tiny-field oracle tests for the two ANF encodings, plus a
reproduction of Trimoska's published n=15, l=5 instance and an end-to-end
WDSat enumeration check.  Run from anywhere:  python3 test_weil.py -v
"""
import os
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import weil  # noqa: E402
import curves  # noqa: E402
import field  # noqa: E402

CACHE = Path(os.environ.get('EC_BASELINE_CACHE', '/tmp/ec-baseline-cache'))
WDSAT_SRC = CACHE / 'vendor/WDSat/src'


def tinyField(n):
    poly, _ = curves.findIrreduciblePoly(n)
    return field.Pb(n, poly)


def oracleRelation(pb, R, xs):
    """Ground truth: some signs give P1 + P2 + P3 = +-R (all x_i must lift)."""
    v = weil.verifyProjected(pb, R, xs)
    return v['lifts'] and v['relation']


def s4Numeric(pb, xs, xR):
    """S4 on the four abscissae, numerically (frozen solver_15 form)."""
    f = pb
    all4 = xs + [xR]
    e1 = 0
    for x in all4:
        e1 ^= x
    e3 = 0
    for i in range(4):
        t = 1
        for j in range(4):
            if i != j:
                t = f.mul(t, all4[j])
        e3 ^= t
    e4 = 1
    for x in all4:
        e4 = f.mul(e4, x)
    sq = lambda v: f.mul(v, v)
    return sq(sq(e3)) ^ f.mul(e4 ^ 1, sq(e3)) ^ f.mul(f.mul(e4, e4 ^ 1), sq(e1)) ^ sq(sq(e1))


def assignmentS4(sysm, pb, d, xs):
    """Assignment of every S'4 variable from the abscissae (aux by definition)."""
    a = {}
    for blk in range(3):
        for k in range(d):
            a[1 + blk * d + k] = (xs[blk] >> k) & 1
    # aux are defined by the first d + (2d-1) + (3d-2) equations: solve them in order
    nDef = d + (2 * d - 1) + (3 * d - 2)
    for i in range(nDef):
        p = sysm.equations[i]
        const, cols = sysm.mono.polyTerms(p)
        auxVar = 3 * d + 1 + i
        rest = const
        for c in cols:
            key = sysm.mono.keys[c]
            if key == (auxVar,):
                continue
            t = 1
            for v in key:
                if not a.get(v, 0):
                    t = 0
                    break
            rest ^= t
        a[auxVar] = rest
    return a


def evalAll(sysm, a, start=0):
    return [weil.evalPoly(sysm.mono, p, a) for p in sysm.equations[start:]]


class S3FormulaTest(unittest.TestCase):
    def test_s4_vanishes_exactly_on_relations_n11(self):
        pb = tinyField(11)
        curve = curves.CurvePb(pb)
        rng = random.Random(1)
        adm = weil.admissibleAbscissae(pb, 11)
        for _ in range(40):
            xs = [rng.choice(adm) for _ in range(3)]
            xR = rng.choice(adm)
            R = curve.pointFromX(xR)
            self.assertEqual(s4Numeric(pb, xs, xR) == 0, oracleRelation(pb, R, xs) or
                             # S4 also vanishes on twist (non-lifting) configurations, but all
                             # four abscissae lift here, so the only other zero is degenerate
                             self._degenerateZero(pb, xs, xR))

    def _degenerateZero(self, pb, xs, xR):
        curve = curves.CurvePb(pb)
        # repeated abscissa with opposite signs cancels: then S4 = S2-type identity
        all4 = xs + [xR]
        for i in range(4):
            for j in range(i + 1, 4):
                if all4[i] == all4[j]:
                    others = [all4[k] for k in range(4) if k not in (i, j)]
                    return others[0] == others[1]
        return False


class S4EncodingTest(unittest.TestCase):
    def test_exhaustive_n11_d3(self):
        pb = tinyField(11)
        curve = curves.CurvePb(pb)
        d = 3
        rng = random.Random(7)
        adm = weil.admissibleAbscissae(pb, d)
        self.assertGreater(len(adm), 2)
        targets = [weil.uniformTarget(pb, rng) for _ in range(3)]
        R, xs, _ = weil.plantedTarget(pb, d, rng)
        targets.append(R)
        for R in targets:
            sysm = weil.buildS4(pb, d, R[0])
            sols = set()
            for x1 in range(1 << d):
                for x2 in range(1 << d):
                    for x3 in range(1 << d):
                        xs = [x1, x2, x3]
                        a = assignmentS4(sysm, pb, d, xs)
                        vals = evalAll(sysm, a)
                        zero = not any(vals)
                        self.assertEqual(zero, s4Numeric(pb, xs, R[0]) == 0,
                                         'ANF disagrees with numeric S4 at %r' % xs)
                        if zero:
                            sols.add(tuple(xs))
            # every projected solution whose abscissae all lift is a genuine relation
            for xs in sols:
                v = weil.verifyProjected(pb, R, list(xs))
                if v['lifts'] and not v['degenerate']:
                    self.assertTrue(v['relation'], xs)

    def test_planted_is_a_solution(self):
        pb = tinyField(13)
        d = 3
        rng = random.Random(3)
        R, xs, _ = weil.plantedTarget(pb, d, rng)
        sysm = weil.buildS4(pb, d, R[0])
        a = assignmentS4(sysm, pb, d, xs)
        self.assertFalse(any(evalAll(sysm, a)))


def rrNumericSolve(pb, R, xs):
    """Solve the RR norm system numerically for (beta, a) given abscissae.
    Returns list of (beta, a, b) solutions."""
    f = pb
    curve = curves.CurvePb(pb)
    r, s = R
    e1 = xs[0] ^ xs[1] ^ xs[2]
    e2 = f.mul(xs[0], xs[1]) ^ f.mul(xs[0], xs[2]) ^ f.mul(xs[1], xs[2])
    e3 = f.mul(f.mul(xs[0], xs[1]), xs[2])
    T = r ^ e1
    if curve.trace(T):
        return []
    out = []
    b0 = curve.halfTrace(T)
    r2 = f.mul(r, r)
    r4 = f.mul(r2, r2)
    rs = r ^ s
    coef = f.mul(rs, rs) ^ 1
    for beta in (0, 1):
        b = b0 ^ beta
        assert f.mul(b, b) ^ b == T
        # E4: r^2 a^2 = r e3 + r^4 + coef b^2
        rhs = f.mul(r, e3) ^ r4 ^ f.mul(coef, f.mul(b, b))
        a2 = f.mul(rhs, f.inv(r2))
        a = f.pow(a2, 1 << (f.m - 1))   # square root
        assert f.mul(a, a) == a2
        E2 = f.mul(a, a) ^ f.mul(a, b) ^ f.mul(r, e1) ^ e2
        E3 = f.mul(r2, b) ^ f.mul(r, f.mul(a, b)) ^ f.mul(rs, f.mul(b, b)) ^ f.mul(r, e2) ^ e3
        if E2 == 0 and E3 == 0:
            out.append((beta, a, b))
    return out


def assignmentRR(sysm, d, xs, beta, a):
    asg = {}
    for blk in range(3):
        for k in range(d):
            asg[1 + blk * d + k] = (xs[blk] >> k) & 1
    asg[sysm.layout['beta']] = beta
    a0 = sysm.layout['a'][0]
    for k in range(sysm.descent.m):
        asg[a0 + k] = (a >> k) & 1
    return asg


class RREncodingTest(unittest.TestCase):
    def test_exhaustive_n11_d3_matches_numeric_and_oracle(self):
        pb = tinyField(11)
        d = 3
        rng = random.Random(11)
        targets = [weil.uniformTarget(pb, rng) for _ in range(3)]
        R, xs, _ = weil.plantedTarget(pb, d, rng)
        targets.append(R)
        for R in targets:
            sysm = weil.buildRR(pb, d, R[0], R[1])
            projected = set()
            for x1 in range(1 << d):
                for x2 in range(1 << d):
                    for x3 in range(1 << d):
                        xs = [x1, x2, x3]
                        sols = rrNumericSolve(pb, R, xs)
                        for beta, a, b in sols:
                            asg = assignmentRR(sysm, d, xs, beta, a)
                            self.assertFalse(any(evalAll(sysm, asg)), (xs, beta, a))
                            # the symbolic b must reproduce the numeric b
                            self.assertEqual(sysm.descent.evaluate(sysm.bSymbolic, asg), b)
                            projected.add(tuple(xs))
                        # a random wrong a must violate the system
                        aBad = rng.getrandbits(pb.m)
                        if not any(aBad == s_[1] for s_ in sols):
                            asg = assignmentRR(sysm, d, xs, rng.getrandbits(1), aBad)
                            self.assertTrue(any(evalAll(sysm, asg)))
            # every RR projected solution with distinct abscissae is a genuine, LIFTING relation
            for xs in projected:
                v = weil.verifyProjected(pb, R, list(xs))
                if not v['degenerate']:
                    self.assertTrue(v['lifts'] and v['relation'], (R, xs))

    def test_planted_solves_and_rr_subset_of_s4(self):
        pb = tinyField(13)
        d = 3
        rng = random.Random(5)
        R, xs, _ = weil.plantedTarget(pb, d, rng)
        self.assertTrue(rrNumericSolve(pb, R, xs))
        for x1 in range(1 << d):
            for x2 in range(1 << d):
                for x3 in range(1 << d):
                    if rrNumericSolve(pb, R, [x1, x2, x3]):
                        self.assertEqual(s4Numeric(pb, [x1, x2, x3], R[0]), 0)


class TrimoskaReproductionTest(unittest.TestCase):
    """Rebuild Xn15l5-1-S.anf from its INFO file and compare equation sets."""

    def parse(self, text):
        eqs = []
        for line in text.splitlines()[1:]:
            toks = line.split()
            assert toks[0] == 'x' and toks[-1] == '0'
            toks = toks[1:-1]
            const = 1
            terms = set()
            i = 0
            while i < len(toks):
                if toks[i] == 'T':
                    const = 0
                    i += 1
                elif toks[i].startswith('.'):
                    k = int(toks[i][1:])
                    terms.add(tuple(sorted(int(t) for t in toks[i + 1:i + 1 + k])))
                    i += 1 + k
                else:
                    terms.add((int(toks[i]),))
                    i += 1
            eqs.append((const, frozenset(terms)))
        return eqs

    @unittest.skipUnless((CACHE / 'instances/Xn15l5-1-S.anf').exists(), 'pinned corpus not fetched')
    def test_same_equations_as_published(self):
        info = (CACHE / 'instances/INFOn15l5-1-S.dimacs').read_text().splitlines()
        n, l = map(int, info[0].split())
        little = lambda s: int(s[::-1], 2)
        poly, xR = little(info[1]), little(info[2])
        pb = field.Pb(n, poly)
        self.assertTrue(pb.isIrreducible())
        mine = weil.buildS4(pb, l, xR)
        ours = self.parse(mine.anfText())
        theirs = self.parse((CACHE / 'instances/Xn15l5-1-S.anf').read_text())
        self.assertEqual(len(ours), len(theirs))
        # aux-definition equations must match line by line (same layout)
        nDef = l + (2 * l - 1) + (3 * l - 2)
        self.assertEqual(ours[:nDef], theirs[:nDef])
        # the S4 rows are the same set of equations (row order is a basis choice)
        self.assertEqual(set(ours[nDef:]), set(theirs[nDef:]))


class WDSatEndToEndTest(unittest.TestCase):
    @unittest.skipUnless(WDSAT_SRC.exists(), 'pinned WDSat source not fetched')
    def test_find_all_enumeration_matches_oracle(self):
        pb = tinyField(11)
        d = 4
        rng = random.Random(21)
        R, planted, _ = weil.plantedTarget(pb, d, rng)
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            for name, sysm in (('s4', weil.buildS4(pb, d, R[0])), ('rr', weil.buildRR(pb, d, R[0], R[1]))):
                cfg, _ = weil.wdsatConfig(sysm.stats(), findAll=True)
                exe = buildWdsat(cfg, tmp / name)
                inp = tmp / (name + '.anf')
                inp.write_text(sysm.anfText())
                out = subprocess.run([str(exe), '-i', str(inp), '-n', str(pb.m), '-l', str(d), '-m', '3', '-x'],
                                     capture_output=True, text=True, timeout=120).stdout
                found = set()
                for line in out.splitlines():
                    if len(line) == sysm.nvars and set(line) <= {'0', '1'}:
                        found.add(tuple(weil.decodeBlocks(line, d)))
                self.assertIn('UNSAT', out)
                # oracle: exhaustive over the system's own semantics
                expect = set()
                for x1 in range(1 << d):
                    for x2 in range(1 << d):
                        for x3 in range(1 << d):
                            xs = [x1, x2, x3]
                            if name == 's4':
                                if s4Numeric(pb, xs, R[0]) == 0:
                                    expect.add(tuple(xs))
                            elif rrNumericSolve(pb, R, xs):
                                expect.add(tuple(xs))
                self.assertEqual(found, expect, name)
                self.assertIn(tuple(planted), found)


def buildWdsat(configText, dest):
    dest.mkdir(parents=True, exist_ok=True)
    src = dest / 'src'
    src.mkdir(exist_ok=True)
    for p in WDSAT_SRC.iterdir():
        if p.suffix in ('.c', '.h'):
            (src / p.name).write_bytes(p.read_bytes())
    (src / 'config.h').write_text(configText)
    exe = dest / 'wdsat_solver'
    cmd = ['gcc', '-O3', '-w'] + sorted(str(p) for p in src.glob('*.c')) + ['-lm', '-o', str(exe)]
    subprocess.run(cmd, check=True, capture_output=True)
    return exe


if __name__ == '__main__':
    unittest.main()
