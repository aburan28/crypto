#!/usr/bin/env python3
# Map the generated GF(2^131) circuits onto 6-input lookup tables, the way an
# FPGA synthesis tool would, and report what an ECC2K-130 iteration costs as
# fabric rather than as instructions.
#
#   python3 lutmap.py [--json result.json]
#
# The GPU ceiling in ../../THROUGHPUT-CEILING.md is an instruction-issue bound:
# the walk's bit operations share one 64-lane integer pipe.  On an FPGA the same
# bit operations are laid down in space, so the question is not how many issue
# slots they need but how many LUTs they occupy and how deep they stack.  This
# measures both, from the same straight-line programs the CUDA client is built
# from, so the two platforms are costed from one circuit.
#
# Every multiplier is verified against field.Onb before it is counted, on 64
# random inputs at a time, so the numbers below describe circuits that compute
# the right answer.
#
# The mapper is a priority-cut mapper: enumerate up to C six-input cuts per
# node, keep the best by depth (or by area flow), then choose a cover.  That is
# the standard FlowMap/priority-cut construction and it lands within a few
# percent of a hand count (see the cross-check in ../../FPGA-CEILING.md).  It
# does not model the LUT_6_2 dual-output packing a real tool uses, so its counts
# are an upper bound -- published hand-packed multipliers are ~40% smaller.
#
# No type hints, camelCase identifiers, no itertools (project convention).

import argparse
import json
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', 'codegen'))

import build
import field
import ir

M = 131
N = 2 * M + 1
INF = 1 << 30


def verifyMul(prog, roots, m, onb, rng):
    """Check the circuit against the independent field model, 64 lanes at once."""
    xs = []
    ys = []
    for _ in range(64):
        xs.append(onb.randomElement(rng))
        ys.append(onb.randomElement(rng))
    vals = {}
    for i in range(m):
        wa = 0
        wb = 0
        for lane in range(64):
            if (onb.toCoords(xs[lane]) >> i) & 1:
                wa |= 1 << lane
            if (onb.toCoords(ys[lane]) >> i) & 1:
                wb |= 1 << lane
        vals[('a', i)] = wa
        vals[('b', i)] = wb
    out = prog.evaluate(vals, roots)
    for lane in range(64):
        got = 0
        for i in range(m):
            if (out[i] >> lane) & 1:
                got |= 1 << i
        if got != onb.toCoords(onb.mul(xs[lane], ys[lane])):
            raise RuntimeError('multiplier mismatch on lane %d' % lane)


def fold(t):
    """Index of gamma_t after folding by gamma_t = gamma_{-t} = gamma_{t mod n}."""
    t %= N
    return t if t <= M else N - t


def mulConversion(cutoff):
    """Bernstein-Lange: onb -> optimal polynomial basis, Karatsuba, back."""
    p = ir.Prog()
    a = []
    b = []
    for i in range(M):
        a.append(p.addInput('a', i))
    for i in range(M):
        b.append(p.addInput('b', i))
    pa = build.multPrepIr(p, a, M)
    pb = build.multPrepIr(p, b, M)
    hh = build.polyMulIr(p, pa, pb, cutoff)
    return p, build.toOnbIr(p, hh, M)


def mulDirectOnb():
    """Quadratic normal-basis multiplier, the shallow alternative.

    gamma_i gamma_j = gamma_{i+j} + gamma_{i-j}, so the coefficient of gamma_k
    is sum_i a_i (b_{fold(k-i)} + b_{fold(k+i)}).  Three times the gates of the
    conversion form and a fifth of the depth, which is the trade an FPGA wants
    and a GPU does not.
    """
    p = ir.Prog()
    a = [None]
    b = [None]
    for i in range(M):
        a.append(p.addInput('a', i))
    for i in range(M):
        b.append(p.addInput('b', i))
    out = []
    for k in range(1, M + 1):
        terms = []
        for i in range(1, M + 1):
            terms.append(p.andOp(a[i], p.xor(b[fold(k - i)], b[fold(k + i)])))
        out.append(p.xorList(terms))
    return p, out


def gateStats(prog, roots):
    """Two-input gate count and gate depth of the live cone."""
    rc = prog.refCounts(roots)
    depth = [0] * len(prog.ops)
    gates = 0
    for i in range(len(prog.ops)):
        node = prog.ops[i]
        if node is None or rc[i] == 0:
            continue
        gates += 1
        best = 0
        for x in node[1]:
            if depth[x] > best:
                best = depth[x]
        depth[i] = 1 + best
    worst = 0
    for r in roots:
        if r is not None and depth[r] > worst:
            worst = depth[r]
    return gates, worst


def lutMap(prog, roots, k=6, cuts=10, areaOnly=False):
    """Priority-cut k-LUT mapping.  Returns (cover, levelOf, depth)."""
    ops = prog.ops
    n = len(ops)
    rc = prog.refCounts(roots)
    isIn = []
    for i in range(n):
        isIn.append(prog.inputRef[i] is not None)
    arrival = [0] * n
    areaFlow = [0.0] * n
    usable = [None] * n          # cuts offered to fanouts, trivial cut first
    impl = [None] * n            # cuts that could implement this node
    for i in range(n):
        if isIn[i]:
            usable[i] = [frozenset((i,))]
            continue
        if ops[i] is None or rc[i] == 0:
            continue
        args = ops[i][1]
        cand = set()
        if len(args) == 1:
            for l in usable[args[0]]:
                cand.add(l)
        else:
            for l1 in usable[args[0]]:
                for l2 in usable[args[1]]:
                    u = l1 | l2
                    if len(u) <= k:
                        cand.add(u)
        scored = []
        for l in cand:
            d = 0
            f = 1.0
            for x in l:
                if arrival[x] > d:
                    d = arrival[x]
                f += areaFlow[x] / max(1, rc[x])
            d += 1
            if areaOnly:
                scored.append((f, d, len(l), l))
            else:
                scored.append((d, f, len(l), l))
        scored.sort(key=lambda t: (t[0], t[1], t[2]))
        keep = scored[:cuts]
        impl[i] = []
        for t in keep:
            if areaOnly:
                impl[i].append((t[3], t[1], t[0]))
            else:
                impl[i].append((t[3], t[0], t[1]))
        best = impl[i][0][1]
        flow = impl[i][0][2]
        for c in impl[i]:
            if c[1] < best:
                best = c[1]
        for c in impl[i]:
            if c[1] == best and c[2] < flow:
                flow = c[2]
        arrival[i] = best
        areaFlow[i] = impl[i][0][2] if areaOnly else flow
        usable[i] = [frozenset((i,))]
        for c in impl[i]:
            usable[i].append(c[0])
    live = []
    for r in roots:
        if r is not None and not isIn[r]:
            live.append(r)
    depth = 0
    for r in live:
        if arrival[r] > depth:
            depth = arrival[r]
    need = [False] * n
    req = [INF] * n
    for r in live:
        need[r] = True
        req[r] = INF if areaOnly else depth
    cover = {}
    for i in range(n - 1, -1, -1):
        if not need[i] or isIn[i]:
            continue
        pick = None
        for c in impl[i]:
            if c[1] > req[i]:
                continue
            if pick is None or (c[2], len(c[0])) < (pick[2], len(pick[0])):
                pick = c
        cover[i] = pick[0]
        for x in pick[0]:
            if not isIn[x]:
                need[x] = True
                if req[i] - 1 < req[x]:
                    req[x] = req[i] - 1
    levelOf = {}
    keys = sorted(cover)
    for i in keys:
        best = 0
        for x in cover[i]:
            if levelOf.get(x, 0) > best:
                best = levelOf.get(x, 0)
        levelOf[i] = 1 + best
    mapped = 0
    for r in live:
        if levelOf[r] > mapped:
            mapped = levelOf[r]
    return cover, levelOf, mapped


def pipelineFfs(prog, cover, levelOf, roots, perStage):
    """Flip-flops needed if a register bank is inserted every perStage levels."""
    lastUse = {}
    for i in cover:
        for x in cover[i]:
            if levelOf[i] > lastUse.get(x, 0):
                lastUse[x] = levelOf[i]
    depth = 0
    for i in levelOf:
        if levelOf[i] > depth:
            depth = levelOf[i]
    for r in roots:
        if r is not None:
            lastUse[r] = depth + 1
    ffs = 0
    stages = (depth + perStage - 1) // perStage
    for s in range(1, stages):
        edge = s * perStage
        for x in lastUse:
            born = 0 if prog.inputRef[x] is not None else levelOf.get(x, 0)
            if born <= edge < lastUse[x]:
                ffs += 1
    return stages, ffs


def measure(name, prog, roots):
    gates, gateDepth = gateStats(prog, roots)
    cover, levelOf, depth = lutMap(prog, roots)
    coverA, levelA, depthA = lutMap(prog, roots, areaOnly=True)
    row = {
        'circuit': name,
        'gates': gates,
        'gateDepth': gateDepth,
        'lutDepthOpt': len(cover),
        'lutDepthOptLevels': depth,
        'lutAreaOpt': len(coverA),
        'lutAreaOptLevels': depthA,
        'pipelining': {},
    }
    for perStage in (1, 2, 3):
        stages, ffs = pipelineFfs(prog, cover, levelOf, roots, perStage)
        row['pipelining']['levelsPerStage%d' % perStage] = {
            'stages': stages, 'flipFlops': ffs}
    print('%-30s gates %6d  depth %3d | LUT6 %6d @ %2d levels | area-opt %6d @ %2d'
          % (name, gates, gateDepth, row['lutDepthOpt'], depth,
             row['lutAreaOpt'], depthA))
    sys.stdout.flush()
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', default='')
    args = ap.parse_args()

    rng = random.Random(7)
    onb = field.Onb(M)
    rows = []

    prog, roots = mulDirectOnb()
    verifyMul(prog, roots, M, onb, rng)
    rows.append(measure('multiply, direct ONB', prog, roots))

    for cutoff in (131, 66, 33, 17, 9):
        prog, roots = mulConversion(cutoff)
        verifyMul(prog, roots, M, onb, rng)
        rows.append(measure('multiply, conversion cut %d' % cutoff, prog, roots))

    prog = ir.Prog()
    x = []
    for i in range(M):
        x.append(prog.addInput('x', i))
    rows.append(measure('hamming weight', prog, build.hammingIr(prog, x, 8)))

    out = {
        'field': 'GF(2^%d), permuted type-II ONB' % M,
        'mapper': {'k': 6, 'priorityCuts': 10,
                   'note': 'no LUT_6_2 dual-output packing; counts are an upper bound'},
        'verified': 'every multiplier checked against field.Onb on 64 random inputs',
        'circuits': rows,
    }
    if args.json:
        fh = open(args.json, 'w')
        json.dump(out, fh, indent=2)
        fh.write('\n')
        fh.close()
        print('wrote ' + args.json)


main()
