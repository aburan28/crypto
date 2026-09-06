#!/usr/bin/env python3
"""Verify a candidate solution against a *shipped* XOR-CNF instance.

Companion to audit_labels.py.  Where that script decides satisfiability
from our own reimplementation of f3, this one checks a witness against
the benchmark file itself, so a bug in our algebra cannot manufacture a
false result.

The encoding is definitional: SAT variables past the free x-coordinate
bits are all *defined* (monomial variables as ANDs of their factors,
e-variables as XORs of monomials).  So assigning the 3*l free bits and
running unit propagation to fixpoint determines every other variable;
we then check that no clause is violated and nothing is left unassigned.

The 3*l lowest-numbered variables are the free bits, in the order
x_{1,0..l-1}, x_{2,0..l-1}, x_{3,0..l-1} -- the upstream variable
compaction (compute_offsets_X) preserves monomial rank order, and the
degree-1 x-monomials rank below every other monomial.

Usage:
    python3 verify_witness.py <instance.dimacs> <l> <X1,X2,X3>

Example (the mislabelled instance found by audit_labels.py):
    python3 verify_witness.py benchmarks/Xn19l6-19-U.dimacs 6 20,34,41
"""

import sys
# Verify a witness against the SHIPPED XOR-CNF instance by unit propagation.
path, l, xs = sys.argv[1], int(sys.argv[2]), [int(t) for t in sys.argv[3].split(',')]

cnf, xors, maxv = [], [], 0
for line in open(path):
    t = line.split()
    if not t or t[0] == 'p': continue
    if t[0] == 'x':
        lits = [int(v) for v in t[1:] if v != '0']
        xors.append(lits); maxv = max(maxv, max(abs(v) for v in lits))
    else:
        lits = [int(v) for v in t if v != '0']
        cnf.append(lits);  maxv = max(maxv, max(abs(v) for v in lits))

val = [None] * (maxv + 1)
# the 3*l lowest-numbered variables are the free x-coordinate bits, in order
# x_{1,0..l-1}, x_{2,0..l-1}, x_{3,0..l-1}
for i, X in enumerate(xs):
    for d in range(l):
        val[i * l + d + 1] = (X >> d) & 1

def lit(v): return None if val[abs(v)] is None else (val[abs(v)] if v > 0 else 1 - val[abs(v)])

changed = True
while changed:
    changed = False
    for c in cnf:                                  # AND-definition clauses
        vs = [lit(v) for v in c]
        if any(x == 1 for x in vs): continue
        unk = [i for i, x in enumerate(vs) if x is None]
        if len(unk) == 1:
            v = c[unk[0]]; val[abs(v)] = 1 if v > 0 else 0; changed = True
    for c in xors:                                 # XOR clauses
        vs = [lit(v) for v in c]
        unk = [i for i, x in enumerate(vs) if x is None]
        if len(unk) == 1:
            rhs = 1 ^ (sum(x for i, x in enumerate(vs) if i != unk[0]) & 1)
            v = c[unk[0]]; val[abs(v)] = rhs if v > 0 else 1 - rhs; changed = True

unassigned = [v for v in range(1, maxv + 1) if val[v] is None]
bad_cnf = [c for c in cnf if not any(lit(v) == 1 for v in c)]
bad_xor = [c for c in xors if (sum(lit(v) for v in c) & 1) != 1]
print('file            :', path.split('/')[-1])
print('variables       :', maxv, ' clauses: %d CNF + %d XOR' % (len(cnf), len(xors)))
print('unassigned after propagation :', len(unassigned))
print('violated CNF clauses         :', len(bad_cnf))
print('violated XOR clauses         :', len(bad_xor))
print('VERDICT         :', 'SATISFIED by this assignment' if not (bad_cnf or bad_xor or unassigned)
      else 'not satisfied')
