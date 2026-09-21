#!/usr/bin/env python3
"""Reconcile this client's permuted type-II ONB with the basis a cairn job pins.

    python3 cairn_basis.py                       # check against the pinned gamma
    python3 cairn_basis.py --gamma <hex>         # ...or another one
    python3 cairn_basis.py --canon <hex,hex,hex> # name the orbit of a corpus record

cairn's orbit piecework objective (aburan28/cairn, docs/design/orbit-piecework.md)
names an orbit by the least cyclic rotation of the abscissa's coordinates in a
normal basis it pins by a generator gamma, indexed so that coordinate i is the
coefficient of gamma^(2^i) and sigma is a rotation.  This client uses a permuted
type-II ONB in which sigma is SQ_PERM instead.  Same basis set, different index
order -- which is a silent disagreement, because both sides produce well-formed
131-bit strings and simply differ about which orbit a point is on.

What this settles: gamma is one of our own basis elements, the two orders differ
by exactly one bijection, and cairn's name for an orbit is recoverable from the
canon[3] a DpFileRecord already carries.  So the 32-byte corpus format does not
have to change to name an orbit; only the witness is new.  See CAIRN-WITNESS.md.

Derive the index, never hardcode it.  A hardcoded index survives a change to the
generator's basis and keeps producing valid-looking names for the wrong orbits,
which nothing downstream would catch.

No type hints, camelCase identifiers, no itertools (project convention).
"""

import argparse
import os
import re

M = 131
PB_TAPS = (13, 2, 1)
POLY = (1 << M) | (1 << 13) | (1 << 2) | (1 << 1) | 1
MASK = (1 << M) - 1

# cairn's pinned generator, examples/certicom-ecdlp/jobs/ecc2k130.json.
CAIRN_GAMMA = 0x85425218439633925EB07B25668EEE8F


def loadGenerated(path):
    """Z_TO_ONB, SQ_PERM and the instance points out of the generated header."""
    src = open(path).read()
    match = re.search(r"Z_TO_ONB\[131\]\[3\] = \{(.*?)\n\};", src, re.S)
    words = [int(v, 16) for v in re.findall(r"0x([0-9a-fA-F]+)ull", match.group(1))]
    if len(words) != 3 * M:
        raise ValueError("Z_TO_ONB is not %d x 3" % M)
    zToOnb = []
    for i in range(M):
        zToOnb.append(words[3 * i] | (words[3 * i + 1] << 64) | (words[3 * i + 2] << 128))
    perm = re.search(r"SQ_PERM\[131\] = \{(.*?)\n\};", src, re.S).group(1)
    sqPerm = [int(v) for v in perm.replace("\n", "").split(",") if v.strip()]
    if len(sqPerm) != M:
        raise ValueError("SQ_PERM is not %d long" % M)
    points = {}
    for name in ("PX", "PY", "QX", "QY"):
        body = re.search(name + r"\[3\] = \{(.*?)\};", src, re.S).group(1)
        limbs = [int(v, 16) for v in re.findall(r"0x([0-9a-fA-F]+)ull", body)]
        points[name] = limbs[0] | (limbs[1] << 64) | (limbs[2] << 128)
    return zToOnb, sqPerm, points


def polySqr(a):
    r = 0
    i = 0
    while a:
        if a & 1:
            r ^= 1 << (2 * i)
        a >>= 1
        i += 1
    while r > MASK:
        hi = r >> M
        r &= MASK
        for tap in PB_TAPS:
            r ^= hi << tap
        r ^= hi
    return r


def normalRows(gamma):
    """Rows whose parity against a polynomial-basis element gives its cairn
    normal-basis coordinates: the inverse of the matrix whose columns are the
    conjugates gamma^(2^i)."""
    columns = []
    cur = gamma
    for _ in range(M):
        columns.append(cur)
        cur = polySqr(cur)
    left = [0] * M
    right = [0] * M
    for i in range(M):
        row = 0
        for j in range(M):
            if (columns[j] >> i) & 1:
                row |= 1 << j
        left[i] = row
        right[i] = 1 << i
    at = 0
    for col in range(M):
        pivot = -1
        for r in range(at, M):
            if (left[r] >> col) & 1:
                pivot = r
                break
        if pivot < 0:
            return None
        left[at], left[pivot] = left[pivot], left[at]
        right[at], right[pivot] = right[pivot], right[at]
        for r in range(M):
            if r != at and ((left[r] >> col) & 1):
                left[r] ^= left[at]
                right[r] ^= right[at]
        at += 1
    out = [0] * M
    for r in range(M):
        out[(left[r] & -left[r]).bit_length() - 1] = right[r]
    return out


def toCairnCoords(value, rows):
    out = 0
    for i in range(M):
        if bin(value & rows[i]).count("1") & 1:
            out |= 1 << i
    return out


def leastRotation(coords):
    best = coords
    cur = coords
    for _ in range(M - 1):
        cur = ((cur << 1) | (cur >> (M - 1))) & MASK
        if cur < best:
            best = cur
    return best


def clientToCairn(value, perm):
    out = 0
    for i in range(M):
        if (value >> perm[i]) & 1:
            out |= 1 << i
    return out


def challengePoints(path):
    """P and Q in POLYNOMIAL basis, from the generator's own constants.

    The cross-check in main() needs a route to cairn's coordinates that shares
    no machinery with the permutation: these go through gamma's conjugate
    matrix, the generated header's ONB values go through SQ_PERM, and the two
    must land on the same bits.  Comparing the permutation against a conversion
    built out of the permutation would check nothing.
    """
    src = open(path).read()
    out = {}
    for name in ("PX", "PY", "QX", "QY"):
        out[name] = int(re.search(r"CHALLENGE_" + name + r"\s*=\s*0x([0-9A-Fa-f]+)",
                                  src).group(1), 16)
    return out


def clientIndexOfGamma(gamma, zToOnb):
    """Which of our ONB coordinates gamma is. Must be exactly one."""
    value = 0
    for i in range(M):
        if (gamma >> i) & 1:
            value ^= zToOnb[i]
    if bin(value).count("1") != 1:
        raise ValueError("gamma is not one of this client's basis elements "
                         "(popcount %d, expected 1)" % bin(value).count("1"))
    return value.bit_length() - 1


def buildPerm(gamma, zToOnb, sqPerm):
    """perm[i] is the client coordinate holding cairn's coordinate i."""
    index = clientIndexOfGamma(gamma, zToOnb)
    perm = [index]
    for _ in range(M - 1):
        perm.append(sqPerm[perm[-1]])
    if sorted(perm) != list(range(M)):
        raise ValueError("the derived index map is not a bijection")
    return index, perm


def clientCanonical(value, sqPerm):
    """The client's own orbit representative: min over the sigma-orbit."""
    best = value
    cur = value
    for _ in range(M - 1):
        nxt = 0
        for i in range(M):
            if (cur >> i) & 1:
                nxt |= 1 << sqPerm[i]
        cur = nxt
        if cur < best:
            best = cur
    return best


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--generated", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "generated", "eccF131.h"))
    parser.add_argument("--gamma", default=None,
                        help="normal-basis generator in polynomial-basis hex")
    parser.add_argument("--gen-py", dest="genPy", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "codegen", "gen.py"))
    parser.add_argument("--canon", default=None,
                        help="a DpFileRecord canon as three comma-separated hex limbs")
    args = parser.parse_args()

    zToOnb, sqPerm, points = loadGenerated(args.generated)
    challenge = challengePoints(args.genPy)
    gamma = int(args.gamma, 16) if args.gamma else CAIRN_GAMMA
    try:
        index, perm = buildPerm(gamma, zToOnb, sqPerm)
    except ValueError as why:
        # The failure this tool exists to catch, so report it rather than
        # unwinding: a gamma outside our basis set is a job we cannot serve.
        print("cannot reconcile: %s" % why)
        return 1
    rows = normalRows(gamma)
    if rows is None:
        raise SystemExit("gamma is not a normal element")

    print("gamma                 %033x" % gamma)
    print("client basis index    %d  (derived, never hardcoded)" % index)
    print("index map             bijection on 0..%d: yes" % (M - 1))

    ok = True
    for name in ("PX", "PY", "QX", "QY"):
        viaPerm = clientToCairn(points[name], perm)
        viaBasis = toCairnCoords(challenge[name], rows)
        agree = viaPerm == viaBasis
        ok = ok and agree
        print("%-4s weight %3d       permutation agrees with the basis conversion: %s"
              % (name, bin(viaPerm).count("1"), "yes" if agree else "NO"))

    canonPx = clientCanonical(points["PX"], sqPerm)
    fromCanon = leastRotation(clientToCairn(canonPx, perm))
    fromPoint = leastRotation(clientToCairn(points["PX"], perm))
    print("our canon and cairn's name are different members of one orbit: %s"
          % ("yes" if clientToCairn(canonPx, perm) != fromPoint else "no"))
    print("cairn's name is recoverable from canon[3] alone: %s"
          % ("yes" if fromCanon == fromPoint else "NO"))
    ok = ok and fromCanon == fromPoint

    if args.canon:
        limbs = [int(v, 16) for v in args.canon.split(",")]
        value = limbs[0] | (limbs[1] << 64) | (limbs[2] << 128)
        print("orbit                 %x" % leastRotation(clientToCairn(value, perm)))

    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
