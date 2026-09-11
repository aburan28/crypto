#!/usr/bin/env python3
# Read the emitted Verilog back from disk and simulate it against field.Onb.
#
#   python3 testverilog.py [--rtl ../generated/rtl] [--m 131] [--vectors 8]
#
# genverilog.py emits RTL from the same straight-line IR gen.py already checks,
# so what is left to go wrong is the *emission*: a mis-indexed coordinate, an
# operand read from the wrong pipeline stage, a register chain one short, an
# output tied to the wrong wire.  None of those show up in the IR.
#
# So this reads the .v files as text, parses them with a small evaluator that
# knows nothing about the generator, and runs them as synchronous circuits:
# drive one input set per clock, advance the registers, and compare the output
# stream against field.Onb.  It also recovers the latency and the initiation
# interval from the simulated waveform rather than trusting the module header,
# which is the pair of claims an interleaved datapath lives on.
#
# The last check is the one that matters most: ecc_pre and ecc_post are driven
# back to back on a real curve point, with the inversion done in Python where
# the batch inverter would do it in hardware, and the result compared against
# sigma^j(R) + R computed by curves.Curve.  That is one whole walk step of the
# datapath, checked against the model the CPU and CUDA clients are checked
# against.
#
# Field elements have a redundant coordinate representation (the all-ones vector
# is zero), so element comparisons go through fromCoords/toCoords rather than
# comparing raw words.  Weights do not: HW is taken of the stored coordinates,
# exactly as the client does, so those are compared as integers.
#
# This is not a substitute for a real simulator -- it covers the subset of
# Verilog the generator emits and nothing else.  tb_ecc_mul131.v is there for
# Icarus, Verilator or Vivado when a toolchain exists.
#
# No type hints, camelCase identifiers, no itertools (project convention).

import argparse
import os
import random
import re

import curves
import field

WIRE = re.compile(r'^\s*wire\s+(w\d+)\s*=\s*(.+);\s*$')
REGA = re.compile(r'^\s*(r\d+_\d+)\s*<=\s*(\S+);\s*$')
OUTA = re.compile(r"^\s*assign\s+(\w+)\[(\d+)\]\s*=\s*(\S+);\s*$")
LAT = re.compile(r'latency (\d+) clocks')
PERM = re.compile(r"^\s*assign\s+(\w+)\[(\d+)\]\s*=\s*(\w+)\[(\d+)\];\s*$")
MUX = re.compile(r"^\s*wire\s+\[\d+:0\]\s+(\w+)\s*=\s*sel\[(\d)\]\s*\?\s*(\w+)\s*:\s*(\w+);\s*$")


def parseModule(path):
    """Parse the subset of Verilog genverilog.py emits."""
    assigns = []
    regs = []
    outputs = {}
    latency = None
    fh = open(path)
    for line in fh:
        if latency is None:
            got = LAT.search(line)
            if got:
                latency = int(got.group(1))
        got = WIRE.match(line)
        if got:
            assigns.append((got.group(1), tokenize(got.group(2))))
            continue
        got = REGA.match(line)
        if got:
            regs.append((got.group(1), got.group(2)))
            continue
        got = OUTA.match(line)
        if got:
            outputs.setdefault(got.group(1), {})[int(got.group(2))] = got.group(3)
    fh.close()
    return assigns, regs, outputs, latency


def tokenize(expr):
    return expr.replace('(', ' ( ').replace(')', ' ) ').split()


def evalExpr(tok, vals):
    """Evaluate one emitted right-hand side.  1-bit values only.

    The grammar is the subset genverilog.py emits: ~, &, ^ and | over one-bit
    nets, with Verilog precedence (~ binds tightest, then &, then ^, then |)
    and explicit parentheses around the fused forms.  Parsing it properly
    rather than matching shapes keeps this honest when the IR grows an op --
    the majority gate in the weight tree is what broke shape matching.
    """
    pos = [0]
    out = parseOr(tok, pos, vals)
    if pos[0] != len(tok):
        raise ValueError('unparsed expression: %s' % ' '.join(tok))
    return out


def parsePrimary(tok, pos, vals):
    if pos[0] >= len(tok):
        raise ValueError('expression ended early: %s' % ' '.join(tok))
    head = tok[pos[0]]
    if head == '(':
        pos[0] += 1
        out = parseOr(tok, pos, vals)
        if pos[0] >= len(tok) or tok[pos[0]] != ')':
            raise ValueError('missing ) in %s' % ' '.join(tok))
        pos[0] += 1
        return out
    if head == '~':
        pos[0] += 1
        return 1 - parsePrimary(tok, pos, vals)
    pos[0] += 1
    return lookup(head, vals)


def parseAnd(tok, pos, vals):
    out = parsePrimary(tok, pos, vals)
    while pos[0] < len(tok) and tok[pos[0]] == '&':
        pos[0] += 1
        out = out & parsePrimary(tok, pos, vals)
    return out


def parseXor(tok, pos, vals):
    out = parseAnd(tok, pos, vals)
    while pos[0] < len(tok) and tok[pos[0]] == '^':
        pos[0] += 1
        out = out ^ parseAnd(tok, pos, vals)
    return out


def parseOr(tok, pos, vals):
    out = parseXor(tok, pos, vals)
    while pos[0] < len(tok) and tok[pos[0]] == '|':
        pos[0] += 1
        out = out | parseXor(tok, pos, vals)
    return out


def lookup(name, vals):
    if name == "1'b0":
        return 0
    if name == "1'b1":
        return 1
    return vals[name]


def simulate(assigns, regs, outputs, wants, drive, cycles):
    """Run the circuit for `cycles` clocks.

    `drive(t)` returns the input bits for clock t; `wants` is [(port, width)].
    Returns one dict of port -> word per clock.
    """
    vals = {}
    for name, _ in regs:
        vals[name] = 0
    trace = []
    for t in range(cycles):
        vals.update(drive(t))
        for name, tok in assigns:
            vals[name] = evalExpr(tok, vals)
        sample = {}
        for port, width in wants:
            word = 0
            for k in range(width):
                src = outputs[port].get(k)
                if src is not None and lookup(src, vals):
                    word |= 1 << k
            sample[port] = word
        trace.append(sample)
        nxt = {}
        for dst, src in regs:
            nxt[dst] = lookup(src, vals)
        vals.update(nxt)
    return trace


def findLatency(trace, expected):
    """Recover the pipeline latency from the waveform, and check II = 1."""
    n = len(expected)
    for d in range(len(trace) - n + 1):
        ok = True
        for i in range(n):
            for port in expected[i]:
                if trace[d + i][port] != expected[i][port]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return d
    return None


def driveWords(pairs):
    """pairs is [(portName, width, word)] -> the flat bit dict the sim wants."""
    bits = {}
    for port, width, word in pairs:
        for i in range(width):
            bits['%s[%d]' % (port, i)] = (word >> i) & 1
    return bits


def sameElement(onb, word, element):
    """Compare a coordinate word against a field element, up to representation."""
    return onb.toCoords(onb.fromCoords(word)) == onb.toCoords(element)


def checkMul(rtl, m, onb, rng, count):
    path = os.path.join(rtl, 'ecc_mul%d.v' % m)
    assigns, regs, outputs, claimed = parseModule(path)
    vectors = []
    for _ in range(count):
        u = onb.randomElement(rng)
        v = onb.randomElement(rng)
        vectors.append((onb.toCoords(u), onb.toCoords(v),
                        onb.toCoords(onb.mul(u, v))))

    def drive(t):
        row = vectors[t] if t < len(vectors) else (0, 0, 0)
        return driveWords([('a', m, row[0]), ('b', m, row[1])])

    trace = simulate(assigns, regs, outputs, [('p', m)], drive,
                     count + claimed + 4)
    expected = []
    for row in vectors:
        expected.append({'p': row[2]})
    got = findLatency(trace, expected)
    if got is None:
        raise RuntimeError('ecc_mul%d: product stream never matched field.Onb' % m)
    if got != claimed:
        raise RuntimeError('ecc_mul%d: latency is %d clocks, header claims %d'
                           % (m, got, claimed))
    print('ecc_mul%d.v: PASS -- %d multiplications, II=1, latency %d clk, '
          'checked against field.Onb' % (m, count, got))


def checkHamming(rtl, m, nBits, rng, count):
    path = os.path.join(rtl, 'ecc_hamming%d.v' % m)
    assigns, regs, outputs, claimed = parseModule(path)
    vectors = []
    for _ in range(count):
        vectors.append(rng.getrandbits(m))

    def drive(t):
        word = vectors[t] if t < len(vectors) else 0
        return driveWords([('x', m, word)])

    trace = simulate(assigns, regs, outputs, [('hw', nBits)], drive,
                     count + claimed + 4)
    expected = []
    for v in vectors:
        expected.append({'hw': bin(v).count('1')})
    got = findLatency(trace, expected)
    if got is None:
        raise RuntimeError('ecc_hamming%d: weights never matched popcount' % m)
    print('ecc_hamming%d.v: PASS -- %d weights, II=1, latency %d clk'
          % (m, count, got))


def preReference(onb, m, nBits, xWord, yWord):
    """What ecc_pre must produce: weight, then d and e for sigma^j."""
    weight = bin(xWord).count('1')
    j = 3 + ((weight // 2) % 8)
    xu = onb.fromCoords(xWord)
    yu = onb.fromCoords(yWord)
    d = onb.add(xu, onb.frob(xu, j))
    e = onb.add(yu, onb.frob(yu, j))
    return weight, j, d, e


def checkPre(rtl, m, nBits, onb, rng, count):
    path = os.path.join(rtl, 'ecc_pre%d.v' % m)
    assigns, regs, outputs, claimed = parseModule(path)
    vectors = []
    for _ in range(count):
        xw = onb.toCoords(onb.randomElement(rng))
        yw = onb.toCoords(onb.randomElement(rng))
        vectors.append((xw, yw))

    def drive(t):
        row = vectors[t] if t < len(vectors) else (0, 0)
        return driveWords([('x', m, row[0]), ('y', m, row[1])])

    wants = [('d', m), ('e', m), ('hw', nBits)]
    trace = simulate(assigns, regs, outputs, wants, drive, count + claimed + 4)
    # the weight decides sigma^j, so check it first and the addends through it
    expected = []
    for xw, yw in vectors:
        weight, _, d, e = preReference(onb, m, nBits, xw, yw)
        expected.append({'hw': weight})
    got = findLatency(trace, expected)
    if got is None:
        raise RuntimeError('ecc_pre%d: weights never matched popcount' % m)
    for i in range(count):
        xw, yw = vectors[i]
        weight, j, d, e = preReference(onb, m, nBits, xw, yw)
        sample = trace[got + i]
        if not sameElement(onb, sample['d'], d):
            raise RuntimeError('ecc_pre%d: vector %d denominator wrong (j=%d)'
                               % (m, i, j))
        if not sameElement(onb, sample['e'], e):
            raise RuntimeError('ecc_pre%d: vector %d numerator wrong (j=%d)'
                               % (m, i, j))
    print('ecc_pre%d.v: PASS -- %d steps, II=1, latency %d clk, weight and '
          'sigma^j addends match field.Onb, 0 multiplications'
          % (m, count, got))


def checkPost(rtl, m, onb, curve, rng, count):
    """Drive ecc_post on real curve points and check sigma^j(R) + R."""
    path = os.path.join(rtl, 'ecc_post%d.v' % m)
    assigns, regs, outputs, claimed = parseModule(path)
    vectors = []
    while len(vectors) < count:
        point = curve.pointFromX(onb.randomElement(rng))
        if point is None:
            continue
        xu, yu = point
        xw = onb.toCoords(xu)
        weight = bin(xw).count('1')
        j = 3 + ((weight // 2) % 8)
        other = curve.frob(point, j)
        d = onb.add(xu, other[0])
        if onb.toCoords(d) == 0:
            continue
        e = onb.add(yu, other[1])
        want = curve.add(point, other)
        if want is None:
            continue
        vectors.append((xw, onb.toCoords(yu), onb.toCoords(d), onb.toCoords(e),
                        onb.toCoords(onb.inv(d)), want))

    def drive(t):
        row = vectors[t] if t < len(vectors) else (0, 0, 0, 0, 0, None)
        return driveWords([('x', m, row[0]), ('y', m, row[1]),
                           ('d', m, row[2]), ('e', m, row[3]),
                           ('di', m, row[4])])

    wants = [('x3', m), ('y3', m)]
    trace = simulate(assigns, regs, outputs, wants, drive, count + claimed + 4)
    expected = []
    for row in vectors:
        expected.append({'x3': onb.toCoords(row[5][0])})
    got = findLatency(trace, expected)
    if got is None:
        raise RuntimeError('ecc_post%d: x3 never matched curves.Curve' % m)
    for i in range(count):
        want = vectors[i][5]
        sample = trace[got + i]
        if not sameElement(onb, sample['x3'], want[0]):
            raise RuntimeError('ecc_post%d: vector %d x3 wrong' % (m, i))
        if not sameElement(onb, sample['y3'], want[1]):
            raise RuntimeError('ecc_post%d: vector %d y3 wrong' % (m, i))
    print('ecc_post%d.v: PASS -- %d additions, II=1, latency %d clk, points '
          'match curves.Curve, 2 multiplications' % (m, count, got))
    return got


def checkStep(rtl, m, nBits, onb, curve, rng, count):
    """One whole walk step through both blocks: R -> sigma^j(R) + R."""
    preA, preR, preO, preLat = parseModule(os.path.join(rtl, 'ecc_pre%d.v' % m))
    postA, postR, postO, postLat = parseModule(os.path.join(rtl, 'ecc_post%d.v' % m))
    points = []
    while len(points) < count:
        point = curve.pointFromX(onb.randomElement(rng))
        if point is None:
            continue
        weight = bin(onb.toCoords(point[0])).count('1')
        j = 3 + ((weight // 2) % 8)
        other = curve.frob(point, j)
        if onb.toCoords(onb.add(point[0], other[0])) == 0:
            continue
        want = curve.add(point, other)
        if want is None:
            continue
        points.append((point, want))

    def drivePre(t):
        row = points[t][0] if t < len(points) else (onb.zero(), onb.zero())
        return driveWords([('x', m, onb.toCoords(row[0])),
                           ('y', m, onb.toCoords(row[1]))])

    preTrace = simulate(preA, preR, preO, [('d', m), ('e', m), ('hw', nBits)],
                        drivePre, count + preLat + 4)
    # the inversion is the one part the batch inverter would do in hardware
    stream = []
    for i in range(count):
        sample = preTrace[preLat + i]
        d = onb.fromCoords(sample['d'])
        stream.append((points[i][0], sample['d'], sample['e'],
                       onb.toCoords(onb.inv(d))))

    def drivePost(t):
        if t < len(stream):
            point, dw, ew, diw = stream[t]
            return driveWords([('x', m, onb.toCoords(point[0])),
                               ('y', m, onb.toCoords(point[1])),
                               ('d', m, dw), ('e', m, ew), ('di', m, diw)])
        return driveWords([('x', m, 0), ('y', m, 0), ('d', m, 0), ('e', m, 0),
                           ('di', m, 0)])

    postTrace = simulate(postA, postR, postO, [('x3', m), ('y3', m)], drivePost,
                         count + postLat + 4)
    for i in range(count):
        sample = postTrace[postLat + i]
        want = points[i][1]
        if not sameElement(onb, sample['x3'], want[0]):
            raise RuntimeError('walk step %d: x wrong' % i)
        if not sameElement(onb, sample['y3'], want[1]):
            raise RuntimeError('walk step %d: y wrong' % i)
        if not curve.onCurve((onb.fromCoords(sample['x3']),
                              onb.fromCoords(sample['y3']))):
            raise RuntimeError('walk step %d: result is not on the curve' % i)
    print('walk step: PASS -- %d iterations of sigma^j(R)+R through ecc_pre and '
          'ecc_post, every result on the curve and equal to curves.Curve'
          % count)


def checkSigma(rtl, m, onb, rng, trials):
    """The Frobenius module is pure wiring; check the permutation it wires."""
    path = os.path.join(rtl, 'ecc_sigma%d.v' % m)
    perms = {}
    muxes = []
    fh = open(path)
    for line in fh:
        got = PERM.match(line)
        if got:
            dst, di, src, si = got.group(1), int(got.group(2)), got.group(3), int(got.group(4))
            perms.setdefault(dst, {})[di] = (src, si)
            continue
        got = MUX.match(line)
        if got:
            muxes.append((got.group(1), int(got.group(2)), got.group(3), got.group(4)))
    fh.close()
    if not perms or len(muxes) != 3:
        raise RuntimeError('ecc_sigma%d: expected one permutation set and 3 muxes' % m)

    def valueOf(netName, sel, xBits, cache):
        if netName == 'x':
            return xBits
        if netName in cache:
            return cache[netName]
        for name, bit, onNet, offNet in muxes:
            if name == netName:
                chosen = onNet if (sel >> bit) & 1 else offNet
                out = valueOf(chosen, sel, xBits, cache)
                cache[netName] = out
                return out
        table = perms[netName]
        src = [None] * m
        for di in range(m):
            srcName, si = table[di]
            src[di] = valueOf(srcName, sel, xBits, cache)[si]
        cache[netName] = src
        return src

    for _ in range(trials):
        u = onb.randomElement(rng)
        coords = onb.toCoords(u)
        xBits = []
        for i in range(m):
            xBits.append((coords >> i) & 1)
        for sel in range(8):
            out = valueOf('c2', sel, xBits, {})
            got = 0
            for i in range(m):
                if out[i]:
                    got |= 1 << i
            want = u
            for _ in range(3 + sel):
                want = onb.mul(want, want)
            if got != onb.toCoords(want):
                raise RuntimeError('ecc_sigma%d: sel=%d is not sigma^%d'
                                   % (m, sel, 3 + sel))
    print('ecc_sigma%d.v: PASS -- sigma^(3+sel) for all 8 selects, %d elements, '
          '0 gates' % (m, trials))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--rtl', default=os.path.join('..', 'generated', 'rtl'))
    ap.add_argument('--m', type=int, default=131)
    ap.add_argument('--vectors', type=int, default=8)
    args = ap.parse_args()

    m = args.m
    rng = random.Random(0xF00D0000 + m)
    onb = field.Onb(m)
    curve = curves.Curve(onb)
    nBits = 1
    while (1 << nBits) <= m:
        nBits += 1

    checkMul(args.rtl, m, onb, rng, args.vectors)
    checkHamming(args.rtl, m, nBits, rng, args.vectors)
    checkSigma(args.rtl, m, onb, rng, 4)
    checkPre(args.rtl, m, nBits, onb, rng, args.vectors)
    checkPost(args.rtl, m, onb, curve, rng, args.vectors)
    checkStep(args.rtl, m, nBits, onb, curve, rng, 4)
    print('all emitted modules agree with field.Onb and curves.Curve')


if __name__ == '__main__':
    main()
