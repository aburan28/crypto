"""Offline search over the walk kernel's build space.

What this can and cannot see is the whole design.  ptxas is deterministic and
needs no GPU, so registers, spill traffic, stack frame and the occupancy that
follows from them are all measurable here.  Wall-clock is not.  So this does not
pick a winner -- it produces the few configurations worth measuring, and
`modal run modal_app.py::autotune` measures those instead of a full cross
product.  Six builds on a real card beats a hundred and forty-four.

The tension it exists to map is occupancy against spilling.  Asking ptxas for
more resident blocks makes it use fewer registers per thread, and past a point
it pays for them by spilling to local memory -- which for this kernel is already
the dominant cost.  Where that trade turns is a curve, and the curve is
something an offline tool can draw even though the optimum on it is not.

Nothing here is a substitute for a measurement.  The session that produced it
had a static analysis say a 66-word leaf was better while the only hardware
available said it was worse, for a reason -- register file size -- that does not
apply to the target.  Treat the ranking as a shortlist, never as a result.

    python3 autolab.py --threads 64,128,256 --min-blocks 1,2,3,4 --leaf 0,17,66

Needs a CUDA include tree and ptxas; see the check-cuda target in the Makefile
for the pip route that supplies both without a toolkit install.

No type hints, camelCase identifiers, no itertools (project convention).
"""

import argparse
import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
REGS_PER_SM = 65536
MAX_WARPS_PER_SM = 64


def run(cmd, cwd=None, timeout=1800):
    p = subprocess.Popen(cmd, shell=True, cwd=cwd, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, text=True)
    try:
        out = p.communicate(timeout=timeout)[0]
    except subprocess.TimeoutExpired:
        p.kill()
        return 1, 'timeout'
    return p.returncode, out


def regenerate(leaf, regs, cudaPath):
    arg = ('--leaf %d' % leaf) if leaf else ('--regs %d' % regs)
    rc, out = run('python3 gen.py --out ../generated %s' % arg, cwd=HERE)
    return rc == 0, out


def compile(batch, threads, minBlocks, arch, cudaPath, ptxas, clang):
    ptx = '/tmp/autolab.ptx'
    defs = '-DECC_BATCH=%d -DECC_THREADS=%d -DECC_MINBLOCKS=%d' % (batch, threads, minBlocks)
    rc, out = run('%s -x cuda --cuda-device-only --cuda-gpu-arch=sm_90 '
                  '--cuda-path=%s -nocudalib -Wno-unknown-cuda-version -O3 -std=c++17 '
                  '%s -S -o %s src/ptxspike.cu' % (clang, cudaPath, defs, ptx), cwd=ROOT)
    if rc != 0:
        return None, 'clang: ' + out[-400:]
    # clang 18 will not emit sm_120; the PTX is architecture neutral for this
    # code, so retarget it and let ptxas do the real architecture work
    src = open(ptx).read()
    src = re.sub(r'^\.version 8\.3', '.version 8.7', src, flags=re.M)
    src = re.sub(r'^\.target sm_90', '.target sm_%s' % arch, src, flags=re.M)
    open(ptx, 'w').write(src)
    rc, log = run('%s -arch=sm_%s -O3 -v %s -o /tmp/autolab.cubin' % (ptxas, arch, ptx))
    if rc != 0:
        return None, 'ptxas: ' + log[-400:]
    return (src, log), None


def parseMetrics(src, log, threads):
    m = re.search(r"Compiling entry function '_Z13eccWalkKernel.*?"
                  r"(\d+) bytes stack frame, (\d+) bytes spill stores, (\d+) bytes spill loads"
                  r".*?Used (\d+) registers", log, re.S)
    if not m:
        return None
    stack, spillSt, spillLd, regs = (int(x) for x in m.groups())
    # Split into functions and match the name, not a substring: the file header
    # names the kernel too, and taking the first hit measured a comment block.
    body = None
    moduleInstrs = moduleLocal = 0
    for part in re.split(r"\n(?=\.(?:visible|weak) \.(?:entry|func))", src):
        head = re.match(r"\.(?:visible|weak) \.(?:entry|func)\s*"
                        r"(?:\([^)]*\)\s*)?([A-Za-z0-9_]+)", part)
        if head is None:
            continue
        n = len(re.findall(r"\n\t(?!//)[a-z]", part))
        l = len(re.findall(r"ld\.local", part)) + len(re.findall(r"st\.local", part))
        moduleInstrs += n
        moduleLocal += l
        if head.group(1).startswith('_Z13eccWalkKernel'):
            body = part
    instrs = len(re.findall(r"\n\t(?!//)[a-z]", body)) if body else 0
    localOps = (len(re.findall(r"ld\.local", body)) +
                len(re.findall(r"st\.local", body))) if body else 0
    blocks = REGS_PER_SM // (regs * threads) if regs * threads else 0
    warps = blocks * threads // 32
    return {'registers': regs, 'stack': stack, 'spillStores': spillSt,
            'spillLoads': spillLd, 'spillBytes': spillSt + spillLd,
            'kernelInstrs': instrs, 'kernelLocalOps': localOps,
            'moduleInstrs': moduleInstrs, 'moduleLocalOps': moduleLocal,
            'blocksPerSM': blocks, 'warpsPerSM': warps,
            'occupancy': round(100.0 * warps / MAX_WARPS_PER_SM, 1)}


def dominates(a, b):
    """a dominates b: at least as good everywhere, strictly better somewhere.

    Better means more warps resident, fewer instructions, less traffic."""
    ge = (a['warpsPerSM'] >= b['warpsPerSM'] and
          a['kernelInstrs'] <= b['kernelInstrs'] and
          a['spillBytes'] + a['kernelLocalOps'] * 4 <= b['spillBytes'] + b['kernelLocalOps'] * 4)
    gt = (a['warpsPerSM'] > b['warpsPerSM'] or
          a['kernelInstrs'] < b['kernelInstrs'] or
          a['spillBytes'] + a['kernelLocalOps'] * 4 < b['spillBytes'] + b['kernelLocalOps'] * 4)
    return ge and gt


def paretoFront(rows):
    out = []
    for r in rows:
        beaten = False
        for s in rows:
            if s is not r and dominates(s, r):
                beaten = True
                break
        if not beaten:
            out.append(r)
    return out


def heuristicCost(r):
    """Work a warp must issue, divided by the warps available to hide it.

    A proxy, and named one: it weights a spilled word at four instructions and
    assumes latency hiding scales with occupancy.  Its job is to order a
    shortlist, not to predict a time."""
    work = r['kernelInstrs'] + r['spillBytes'] / 4.0 + r['kernelLocalOps']
    return work / max(1, r['warpsPerSM'])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--batch', default='8,32')
    ap.add_argument('--threads', default='64,128,256')
    ap.add_argument('--min-blocks', default='1,2,3,4')
    ap.add_argument('--leaf', default='0,17',
                    help='generated leaf sizes; 0 means the register-budget choice')
    ap.add_argument('--regs', type=int, default=255)
    ap.add_argument('--arch', default='120')
    ap.add_argument('--cuda-path', default='/tmp/cudaroot')
    ap.add_argument('--ptxas', default='/tmp/cudaroot/bin/ptxas')
    ap.add_argument('--clang', default='clang')
    ap.add_argument('--out', default='autolab.json')
    ap.add_argument('--top', type=int, default=6)
    args = ap.parse_args()

    ints = lambda s: [int(x) for x in s.split(',') if x != '']
    batches, threadList = ints(args.batch), ints(args.threads)
    minBlocks, leaves = ints(args.min_blocks), ints(args.leaf)

    cache = {}
    if os.path.exists(args.out):
        cache = json.load(open(args.out))

    total = len(batches) * len(threadList) * len(minBlocks) * len(leaves)
    print('%d configurations, about %.0f min at 23 s each' % (total, total * 23 / 60.0))
    rows = []
    n = 0
    curLeaf = None
    for leaf in leaves:
        if leaf != curLeaf:
            ok, out = regenerate(leaf, args.regs, args.cuda_path)
            if not ok:
                print('generator failed for leaf %d: %s' % (leaf, out[-300:]))
                continue
            curLeaf = leaf
        for batch in batches:
            for threads in threadList:
                for mb in minBlocks:
                    n += 1
                    key = '%d/%d/%d/%d' % (leaf, batch, threads, mb)
                    if key in cache:
                        rows.append(cache[key])
                        print('[%d/%d] %-16s cached' % (n, total, key), flush=True)
                        continue
                    t0 = time.time()
                    got, err = compile(batch, threads, mb, args.arch,
                                       args.cuda_path, args.ptxas, args.clang)
                    if err:
                        print('[%d/%d] %-16s FAILED %s' % (n, total, key, err[:80]), flush=True)
                        continue
                    met = parseMetrics(got[0], got[1], threads)
                    if met is None:
                        print('[%d/%d] %-16s no walk kernel in the log' % (n, total, key), flush=True)
                        continue
                    met.update({'leaf': leaf, 'batch': batch, 'threads': threads,
                                'minBlocks': mb, 'seconds': round(time.time() - t0, 1)})
                    met['cost'] = round(heuristicCost(met), 1)
                    rows.append(met)
                    cache[key] = met
                    json.dump(cache, open(args.out, 'w'), indent=1)
                    print('[%d/%d] %-16s regs %3d  warps/SM %2d  instrs %6d  '
                          'spill %5dB  local %5d  cost %7.1f'
                          % (n, total, key, met['registers'], met['warpsPerSM'],
                             met['kernelInstrs'], met['spillBytes'],
                             met['kernelLocalOps'], met['cost']), flush=True)

    if not rows:
        print('nothing measured')
        return 1
    front = paretoFront(rows)
    front.sort(key=lambda r: r['cost'])
    print('\n%d configurations, %d on the Pareto front '
          '(occupancy vs instructions vs traffic)' % (len(rows), len(front)))
    print('  %-6s %-6s %-8s %-10s %5s %6s %8s %8s %8s' %
          ('leaf', 'batch', 'threads', 'minBlocks', 'regs', 'warps', 'instrs', 'spillB', 'cost'))
    for r in front[:args.top]:
        print('  %-6d %-6d %-8d %-10d %5d %6d %8d %8d %8.1f'
              % (r['leaf'], r['batch'], r['threads'], r['minBlocks'],
                 r['registers'], r['warpsPerSM'], r['kernelInstrs'],
                 r['spillBytes'], r['cost']))

    short = front[:args.top]
    print('\nMeasure these on a GPU -- the ranking above is a proxy, not a result:')
    print('  modal run modal_app.py::autotune \\')
    print('      --batches %s \\' % ','.join(sorted(set(str(r['batch']) for r in short))))
    print('      --thread-counts %s \\' % ','.join(sorted(set(str(r['threads']) for r in short))))
    print('      --min-blocks-list %s \\' % ','.join(sorted(set(str(r['minBlocks']) for r in short))))
    print('      --leaves %s' % ','.join(sorted(set(str(r['leaf']) for r in short))))
    return 0


if __name__ == '__main__':
    sys.exit(main())
