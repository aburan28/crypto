#!/usr/bin/env python3
"""Audit the SAT/UNSAT labels of the EC-Index-Calculus-Benchmarks corpus.

Upstream: https://github.com/mtrimoska/EC-Index-Calculus-Benchmarks

Each shipped instance is a Weil descent of the symmetrised 4th Semaev
polynomial for the Koblitz curve y^2 + xy = x^3 + x^2 + 1 over GF(2^n),
with the three unknown x-coordinates confined to the l-dimensional
subspace <1, a, ..., a^(l-1)>.  The generator labels an instance "S" when
it planted a decomposition and "U" when it drew a random x_R -- but "U"
is an *expected* outcome, not a certificate.

This script decides every instance exhaustively:

  * "S" instances -- re-evaluate f3 on the planted solution recorded in
    the INFO file and confirm it vanishes.
  * "U" instances -- search every sorted triple (X1 <= X2 <= X3) in the
    subspace.  Sorted triples suffice because f3 is symmetric in
    X1, X2, X3 (it is built from e1, e2, e3).

Result as of this writing: all 30 planted solutions verify, and 1 of the
30 "U" instances (n19l6-19-U) is in fact satisfiable.

Usage:
    python3 audit_labels.py /path/to/EC-Index-Calculus-Benchmarks
"""

import glob, os, itertools, sys

def parse(path):
    L = open(path).read().split('\n')
    n, l = map(int, L[0].split())
    irr = int(L[1][::-1], 2)          # bit k of string = coeff of x^k
    xr  = int(L[2][::-1], 2)
    lab = L[3].strip()
    sol = L[4].strip()
    return n, l, irr, xr, lab, sol

def make_field(n, irr):
    def mul(a, b):
        r = 0
        while b:
            if b & 1: r ^= a
            b >>= 1
            a <<= 1
            if a >> n & 1: a ^= irr
        return r
    return mul

def f3(mul, X1, X2, X3, Xr):
    sq = lambda z: mul(z, z)
    e1 = X1 ^ X2 ^ X3
    e2 = mul(X1, X2) ^ mul(X2, X3) ^ mul(X1, X3)
    e3 = mul(mul(X1, X2), X3)
    e1_2, e2_2, e3_2 = sq(e1), sq(e2), sq(e3)
    e1_4, e2_4, e3_4 = sq(e1_2), sq(e2_2), sq(e3_2)
    e3_3 = mul(e3_2, e3)
    Xr2 = sq(Xr); Xr3 = mul(Xr2, Xr); Xr4 = sq(Xr2)
    return (Xr4 ^ e1_4 ^ e3_4 ^ mul(e2_4, Xr4) ^ mul(e3_3, Xr)
            ^ mul(mul(e3, e2_2), Xr3) ^ mul(mul(e3, e1_2), Xr) ^ mul(e3, Xr3)
            ^ mul(mul(e1_2, e3_2), Xr2) ^ mul(e3_2, Xr4) ^ e3_2 ^ mul(e2_2, Xr2))

base = sys.argv[1]
files = sorted(glob.glob(os.path.join(base, 'benchmarks', 'INFO*.dimacs')))
stats = {}
for p in files:
    n, l, irr, xr, lab, sol = parse(p)
    mul = make_field(n, irr)
    fam = 'n%dl%d' % (n, l)
    st = stats.setdefault(fam, dict(S=0, U=0, S_ok=0, U_sat=0, U_witness=[]))
    st[lab] += 1
    if lab == 'S':                      # verify the planted solution really works
        X1, X2, X3 = [int(s[::-1], 2) for s in sol.split('-')]
        if f3(mul, X1, X2, X3, xr) == 0: st['S_ok'] += 1
        else: print('PLANTED SOLUTION FAILS:', p)
    else:                               # exhaustively decide the "UNSAT" instance
        sat = None
        for X1 in range(1 << l):
            for X2 in range(X1, 1 << l):
                for X3 in range(X2, 1 << l):
                    if f3(mul, X1, X2, X3, xr) == 0:
                        sat = (X1, X2, X3); break
                if sat: break
            if sat: break
        if sat:
            st['U_sat'] += 1
            st['U_witness'].append((os.path.basename(p), sat))

for fam, st in sorted(stats.items()):
    print('%-8s  S=%2d (planted sol verified: %2d)   U=%2d  of which ACTUALLY SATISFIABLE: %d'
          % (fam, st['S'], st['S_ok'], st['U'], st['U_sat']))
    for w in st['U_witness']:
        print('              mislabelled:', w[0], '-> X1,X2,X3 =', w[1])
