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

Two ways to get PTX out.  With a CUDA toolkit installed, `--compiler nvcc` is
the honest one: nvcc emits PTX for the target architecture directly.  Without
one, `--compiler clang` needs only a CUDA include tree and a standalone ptxas
(see the check-cuda target in the Makefile for the pip route that supplies
both), but clang 18 will not emit sm_120, so it compiles for sm_90 and retargets
the PTX header before handing it to ptxas.  That is sound for this kernel, which
uses no architecture-specific instructions, and unsound in general.  Prefer nvcc
where it exists; `modal run modal_app.py::autolab` runs it in a CPU container
that has one.

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
CACHE_SCHEMA = 3


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


def compile(batch, threads, minBlocks, arch, cudaPath, ptxas, clang, compiler):
    ptx = '/tmp/autolab.ptx'
    defs = '-DECC_BATCH=%d -DECC_THREADS=%d -DECC_MINBLOCKS=%d' % (batch, threads, minBlocks)
    if compiler == 'nvcc':
        rc, out = run('%s/bin/nvcc -ptx -arch=compute_%s -O3 -std=c++17 %s '
                      '-o %s src/ptxspike.cu' % (cudaPath, arch, defs, ptx), cwd=ROOT)
        if rc != 0:
            return None, 'nvcc: ' + out[-400:]
        src = open(ptx).read()
    else:
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


# A PTX instruction is a statement and a statement ends in a semicolon, except a
# call, whose argument list both compilers spread over several lines -- so calls
# are counted separately.  Register and parameter declarations start with a dot
# and are not instructions; nor are the bare identifier lines inside a call's
# argument list, which a looser rule counted, inflating a configuration in
# proportion to how many calls it makes -- which is the thing being compared.
INSTR_RE = re.compile(r"\n\t(?!\.)[@a-z][^\n]*;")
CALL_RE = re.compile(r"\bcall\.uni\b")
LOCAL_RE = re.compile(r"(?:ld|st)\.local")


def countPtx(part):
    instrs = len(INSTR_RE.findall(part)) + len(CALL_RE.findall(part))
    return instrs, len(LOCAL_RE.findall(part))


def splitFunctions(src):
    """Map every PTX function name to its text.

    A name appears twice, once as a forward declaration and once as the
    definition; the definition comes second and wins.  Splitting on function
    boundaries and matching the name is what makes this reliable -- a plain
    search for the kernel's name finds it first in the file header comment, and
    then measures the comment."""
    parts = {}
    for part in re.split(r"\n(?=\.(?:visible|weak) \.(?:entry|func))", src):
        head = re.match(r"\.(?:visible|weak) \.(?:entry|func)\s*"
                        r"(?:\([^)]*\)\s*)?([A-Za-z0-9_$]+)", part)
        if head is not None:
            parts[head.group(1)] = part
    return parts


def reachableFrom(parts, entry):
    """Names the entry can reach, following call.uni edges transitively.

    The walk kernel calls out to __noinline__ helpers -- the multiplier above
    all -- so counting only the entry function is blind to exactly the knob this
    tool exists to turn.  Counting the whole module is wrong the other way: it
    charges the walk for the init and reseed kernels, which run once.  The call
    closure is the code a walk step actually issues."""
    seen, stack = set(), [entry]
    while stack:
        n = stack.pop()
        if n in seen or n not in parts:
            continue
        seen.add(n)
        stack.extend(re.findall(r"call\.uni\s*(?:\([^)]*\)\s*,\s*)?([A-Za-z0-9_$]+)\s*,",
                                parts[n]))
    return seen


def parseMetrics(src, log, threads):
    m = re.search(r"Compiling entry function '_Z13eccWalkKernel.*?"
                  r"(\d+) bytes stack frame, (\d+) bytes spill stores, (\d+) bytes spill loads"
                  r".*?Used (\d+) registers", log, re.S)
    if not m:
        return None
    stack, spillSt, spillLd, regs = (int(x) for x in m.groups())
    parts = splitFunctions(src)
    entry = None
    moduleInstrs = moduleLocal = 0
    for name, part in parts.items():
        n, l = countPtx(part)
        moduleInstrs += n
        moduleLocal += l
        if name.startswith('_Z13eccWalkKernel'):
            entry = name
    walkInstrs = walkLocal = kernelInstrs = kernelLocal = 0
    if entry is not None:
        kernelInstrs, kernelLocal = countPtx(parts[entry])
        for name in reachableFrom(parts, entry):
            n, l = countPtx(parts[name])
            walkInstrs += n
            walkLocal += l
    blocks = REGS_PER_SM // (regs * threads) if regs * threads else 0
    warps = blocks * threads // 32
    return {'registers': regs, 'stack': stack, 'spillStores': spillSt,
            'spillLoads': spillLd, 'spillBytes': spillSt + spillLd,
            'walkInstrs': walkInstrs, 'walkLocalOps': walkLocal,
            'kernelInstrs': kernelInstrs, 'kernelLocalOps': kernelLocal,
            'moduleInstrs': moduleInstrs, 'moduleLocalOps': moduleLocal,
            'blocksPerSM': blocks, 'warpsPerSM': warps,
            'occupancy': round(100.0 * warps / MAX_WARPS_PER_SM, 1)}


def dominates(a, b):
    """a dominates b: at least as good everywhere, strictly better somewhere.

    Better means more warps resident, fewer instructions, less traffic."""
    aTraffic = a['spillBytes'] + a['walkLocalOps'] * 4
    bTraffic = b['spillBytes'] + b['walkLocalOps'] * 4
    ge = (a['warpsPerSM'] >= b['warpsPerSM'] and
          a['walkInstrs'] <= b['walkInstrs'] and aTraffic <= bTraffic)
    gt = (a['warpsPerSM'] > b['warpsPerSM'] or
          a['walkInstrs'] < b['walkInstrs'] or aTraffic < bTraffic)
    return ge and gt


def bindingBound(r):
    """Whether __launch_bounds__ actually constrained the register allocator.

    Occupancy here is register limited: an SM has 65536 registers, so the
    resident warps are 2048/registers and the block size does not enter.  The
    bound only bites once minBlocks * threads exceeds 65536/255, about 257
    threads -- below that ptxas takes its 255 and every minBlocks gives an
    identical build.  Saying so keeps eight copies of one configuration off a
    shortlist that exists to be measured on a card."""
    return r['minBlocks'] * r['threads'] > REGS_PER_SM // 255


def dedupe(rows):
    """Collapse configurations that compiled to the same thing."""
    seen = {}
    for r in rows:
        key = (r['leaf'], r['registers'], r['walkInstrs'],
               r['spillBytes'], r['walkLocalOps'], r['warpsPerSM'])
        if key not in seen or (r['minBlocks'], r['threads']) < (seen[key]['minBlocks'],
                                                               seen[key]['threads']):
            seen[key] = r
    return list(seen.values())


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
    shortlist, not to predict a time.

    It counts the call closure rather than the entry function.  An earlier
    version counted only the entry and could not tell a 17-word multiplier leaf
    from the register-budget default, because the multiplier is __noinline__ and
    every instruction that differs between them lives in the callee."""
    work = r['walkInstrs'] + r['spillBytes'] / 4.0 + r['walkLocalOps']
    return work / max(1, r['warpsPerSM'])


SELFTEST_PTX = """//
// Generated by NVIDIA NVVM Compiler
// _Z13eccWalkKernelI7CfgF131jEv10WalkParamsIT0_E is named here, in a comment
//
.version 8.7
.target sm_120
.address_size 64

.weak .func _ZN3mulE;
.visible .entry _Z13eccWalkKernelI7CfgF131jEv10WalkParamsIT0_E(
\t.param .align 8 .b8 p[8]
)
{
\t.reg .b32 %r<4>;
\tld.param.u64 %rd1, [p];
\tst.local.u32 [%rd1], %r1;
\tcall.uni
\t_ZN3mulE,
\t(
\tparam0
\t);
\tret;
}

.weak .func _ZN3mulE(
\t.param .b64 q
)
{
\txor.b32 %r1, %r2, %r3;
\tld.local.u32 %r2, [%rd1];
\tret;
}

.visible .entry _Z13eccInitKernelI7CfgF131jEv10WalkParamsIT0_E()
{
\tmov.u32 %r1, 0;
\tmov.u32 %r2, 0;
\tmov.u32 %r3, 0;
\tret;
}
"""


def selfTest():
    """Guard the extraction bugs that produced plausible wrong numbers.

    All three were silent -- a name matched inside the file header comment
    reported zero instructions, counting only the entry function reported the
    same total for every multiplier leaf, and counting indented lines charged a
    call once per argument.  None of them failed, so none of them was noticed."""
    parts = splitFunctions(SELFTEST_PTX)
    walk = '_Z13eccWalkKernelI7CfgF131jEv10WalkParamsIT0_E'
    init = '_Z13eccInitKernelI7CfgF131jEv10WalkParamsIT0_E'
    assert set(parts) == {walk, init, '_ZN3mulE'}, sorted(parts)
    # the definition must win over the forward declaration
    assert 'xor.b32' in parts['_ZN3mulE']
    assert reachableFrom(parts, walk) == {walk, '_ZN3mulE'}
    # three statements plus one call; the call's argument lines are not counted
    assert countPtx(parts[walk]) == (4, 1), countPtx(parts[walk])
    assert countPtx(parts['_ZN3mulE']) == (3, 1), countPtx(parts['_ZN3mulE'])
    assert countPtx(parts[init]) == (4, 0), countPtx(parts[init])

    log = ("ptxas info    : Compiling entry function "
           "'_Z13eccWalkKernelI7CfgF131jEv10WalkParamsIT0_E' for 'sm_120'\n"
           "ptxas info    : Function properties for "
           "_Z13eccWalkKernelI7CfgF131jEv10WalkParamsIT0_E\n"
           "    16 bytes stack frame, 8 bytes spill stores, 4 bytes spill loads\n"
           "ptxas info    : Used 64 registers\n")
    m = parseMetrics(SELFTEST_PTX, log, 256)
    assert m['registers'] == 64 and m['spillBytes'] == 12, m
    # the entry alone sees 4; the closure sees the multiplier too; the module
    # also charges the init kernel, which runs once and is not the walk's cost
    assert (m['kernelInstrs'], m['walkInstrs'], m['moduleInstrs']) == (4, 7, 11), m
    assert (m['kernelLocalOps'], m['walkLocalOps'], m['moduleLocalOps']) == (1, 2, 2), m
    # 65536/64 registers is 1024 resident threads, which is 32 warps
    assert m['warpsPerSM'] == 32, m
    assert bindingBound({'minBlocks': 4, 'threads': 256})
    assert not bindingBound({'minBlocks': 2, 'threads': 128})
    print('autolab self-test: PTX extraction and occupancy model agree')


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
    ap.add_argument('--compiler', default='auto', choices=('auto', 'nvcc', 'clang'),
                    help='auto picks nvcc when --cuda-path has one')
    ap.add_argument('--out', default='autolab.json')
    ap.add_argument('--top', type=int, default=6)
    ap.add_argument('--self-test', action='store_true',
                    help='check the PTX extraction against a fixture and exit')
    args = ap.parse_args()

    if args.self_test:
        selfTest()
        return 0

    compiler = args.compiler
    if compiler == 'auto':
        compiler = 'nvcc' if os.path.exists(
            os.path.join(args.cuda_path, 'bin', 'nvcc')) else 'clang'
    print('compiling with %s, ptxas targeting sm_%s' % (compiler, args.arch))

    ints = lambda s: [int(x) for x in s.split(',') if x != '']
    batches, threadList = ints(args.batch), ints(args.threads)
    minBlocks, leaves = ints(args.min_blocks), ints(args.leaf)

    # A cached row is only comparable to a fresh build if the same metrics were
    # extracted from it, so key the file by schema and start over when that
    # changes rather than silently mixing two definitions of the cost.
    cache = {}
    if os.path.exists(args.out):
        loaded = json.load(open(args.out))
        if loaded.get('schema') == CACHE_SCHEMA and loaded.get('compiler') == compiler:
            cache = loaded.get('rows', {})
        else:
            print('%s came from a different build of this tool; recompiling' % args.out)

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
                    got, err = compile(batch, threads, mb, args.arch, args.cuda_path,
                                       args.ptxas, args.clang, compiler)
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
                    json.dump({'schema': CACHE_SCHEMA, 'compiler': compiler,
                               'arch': args.arch, 'rows': cache},
                              open(args.out, 'w'), indent=1)
                    print('[%d/%d] %-16s regs %3d  warps/SM %2d  instrs %6d  '
                          'spill %5dB  local %5d  cost %7.1f'
                          % (n, total, key, met['registers'], met['warpsPerSM'],
                             met['walkInstrs'], met['spillBytes'],
                             met['walkLocalOps'], met['cost']), flush=True)

    if not rows:
        print('nothing measured')
        return 1
    unique = dedupe(rows)
    if len(unique) < len(rows):
        loose = len([r for r in rows if not bindingBound(r)])
        print('\n%d of %d configurations compiled to something already seen; %d of '
              'them asked for a launch bound that cannot bind (minBlocks * threads '
              'below %d)' % (len(rows) - len(unique), len(rows), loose,
                             REGS_PER_SM // 255 + 1))
    front = paretoFront(unique)
    front.sort(key=lambda r: r['cost'])
    print('\n%d distinct builds, %d on the Pareto front '
          '(occupancy vs instructions vs traffic)' % (len(unique), len(front)))
    print('  %-6s %-6s %-8s %-10s %5s %6s %8s %8s %8s' %
          ('leaf', 'batch', 'threads', 'minBlocks', 'regs', 'warps', 'instrs', 'spillB', 'cost'))
    for r in front[:args.top]:
        print('  %-6d %-6d %-8d %-10d %5d %6d %8d %8d %8.1f'
              % (r['leaf'], r['batch'], r['threads'], r['minBlocks'],
                 r['registers'], r['warpsPerSM'], r['walkInstrs'],
                 r['spillBytes'], r['cost']))

    # Name the front's rows, not the distinct values in them.  Crossing the
    # values back out measures two or three times as many builds as the front
    # has points, which gives away the whole reason for searching offline first.
    short = front[:args.top]
    plan = ','.join('%d:%d:%d:%d' % (r['leaf'], r['batch'], r['threads'], r['minBlocks'])
                    for r in short)
    print('\nMeasure these %d builds on a GPU -- the ranking above is a proxy, '
          'not a result:' % len(short))
    print('  modal run modal_app.py::autotune --configs %s' % plan)
    return 0


if __name__ == '__main__':
    sys.exit(main())
