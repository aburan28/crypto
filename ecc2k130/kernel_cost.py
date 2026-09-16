#!/usr/bin/env python3
"""Static SASS cost of the packed walk kernel per scalar update.

Compiles packedkernels.cuh's walk() with the audited RTX PRO 6000 preset (or
overrides), dumps sm_120 SASS, finds the step and slot loops from backward
branches, attributes every __noinline__ callee's body to its call sites, and
weights each loop body by its trip count: the two slot loops run ECC_BATCH
times per step, the inversion and step overhead once.  Reports alu / imad /
clmad / mem lane-instructions per update, plus registers and spills from
ptxas -v.

This replaces the per-routine numbers of sass_cost.py for whole-kernel
questions: that tool counts routine bodies in isolation and so counted the
ECC_BIG (__noinline__) products and Frobenius networks zero times inside the
kernel and weighted the paired product as one product.  See
ITERATION-FUNCTION.md section 3.

Static counts include both arms of predicated branches, so they are an upper
bound on the dynamic count; the measured rate bounds the dynamic count from
the other side (ITERATION-FUNCTION.md section 1).

Needs nvcc and cuobjdump (CUDA 13.3+ for clmad) on PATH:
  ./kernel_cost.py                 # audited preset
  ./kernel_cost.py --ops           # plus an opcode histogram of the slot loops
  ./kernel_cost.py --define ECC_PACKED_CLMAD=0
"""
import argparse, os, re, subprocess, sys, tempfile
from collections import Counter

HERE = os.environ.get("ECC2K130", os.path.dirname(os.path.abspath(__file__)))
PRESET = {
    "ECC_BATCH": 16, "ECC_THREADS": 256, "ECC_MINBLOCKS": 2,
    "ECC_PACKED_SINGLE_PRODUCT": 1, "ECC_PACKED_CACHE_DENOM": 1,
    "ECC_PACKED_BY_VALUE": 1, "ECC_PACKED_PERM_SIGMA": 3,
    "ECC_PACKED_POLY_CHAIN": 1, "ECC_PACKED_UNROLL_INV": 1,
    "ECC_PACKED_PAIR_PRODUCTS": 1, "ECC_PACKED_POLY_STATE": 1,
    "ECC_PACKED_DIRECT_REDUCE": 1, "ECC_PACKED_GENERATED_PRODUCT": 1,
    "ECC_PACKED_CLMAD": 1, "ECC_PACKED_WEIGHTED_PREFIX": 2,
    "ECC_PACKED_COMPACT_STATE": 1, "ECC_PACKED_SHARED_SIGMA": 1,
    "ECC_PACKED_STATE_TILE": 256,
}
SPIKE = r'''
#include "../include/curveparams.h"
#include "../include/packedkernels.cuh"
'''
ALU = ("LOP3","LOP","IADD3","IADD","SHF","SHL","SHR","PRMT","POPC","FLO","BREV","SEL",
       "ISETP","IMNMX","MOV","PLOP3","P2R","R2P","VOTE","IABS","LEA","XOR","BFE","BFI",
       "BMSK","FSEL","F2I","I2F","IADD32I","VIADD","VIADDMNMX")
MEM = ("LDG","STG","LDS","STS","LDL","STL","LD","ST","LDC","RED","ATOM","ATOMG","MEMBAR",
       "CCTL","LDSM","ULDC")
CTRL = ("BRA","BSYNC","BSSY","EXIT","CALL","RET","NOP","BAR","SSY","SYNC","JMP","PBK","BRK",
        "WARPSYNC","YIELD","S2R","CS2R","DEPBAR","BMOV","BPT","ERRBAR","USETP","UISETP",
        "UMOV","UIADD3","ULOP3","USHF","R2UR","UPRMT","ULEA","UBMSK")

def classify(op):
    head = op.split(".")[0]
    if head.startswith("CLMAD") or head.startswith("CLMUL"): return "clmad"
    if op.startswith(("IMAD.MOV","IMAD.IADD","IMAD.SHL")): return "alu"
    if head.startswith("IMAD") or head.startswith("IMUL"): return "imad"
    if head in MEM: return "mem"
    if head in CTRL: return "ctrl"
    return "alu"

def build(defs, arch, work, extra):
    src = os.path.join(work, "walk_spike.cu"); open(src,"w").write(SPIKE)
    cubin = os.path.join(work, "walk.cubin")
    cmd = ["nvcc","-O3","-std=c++17","-arch="+arch,"-cubin","-Xptxas","-v",
           "-I",HERE,"-I",os.path.join(HERE,"include"),"-o",cubin,src] + extra
    cmd[1:1] = ["-D%s=%s"%(k,v) for k,v in defs.items()]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=os.path.join(HERE,"src"))
    if r.returncode: sys.exit("nvcc failed:\n"+r.stderr[:4000])
    return cubin, r.stderr

def sass(cubin, kernel):
    r = subprocess.run(["cuobjdump","-sass",cubin], capture_output=True, text=True)
    if r.returncode: sys.exit(r.stderr)
    fn, insts = None, []
    for line in r.stdout.splitlines():
        m = re.match(r"\s*Function : (\S+)", line)
        if m:
            fn = m.group(1); continue
        if fn is None or kernel not in fn: continue
        m = re.match(r"\s*/\*([0-9a-f]{4,})\*/\s+(@!?U?P\d\s+)?([A-Z][A-Z0-9_.]*)(.*)", line)
        if m:
            insts.append((int(m.group(1),16), m.group(3), m.group(4), bool(m.group(2))))
    return insts

def loops(insts):
    """(start_addr, end_addr) for every backward branch."""
    out = []
    for addr, op, rest, pred in insts:
        if op.startswith("BRA"):
            m = re.search(r"0x([0-9a-f]+)", rest)
            if m:
                tgt = int(m.group(1),16)
                if tgt <= addr: out.append((tgt, addr))
    return out

def count(insts, lo, hi):
    c = Counter()
    for addr, op, rest, pred in insts:
        if lo <= addr <= hi: c[classify(op)] += 1; c["n"] += 1
    return c

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arch", default="sm_120")
    ap.add_argument("--define", action="append", default=[])
    ap.add_argument("--kernel", default="4walk")
    ap.add_argument("--dump", help="write SASS listing here")
    ap.add_argument("--extra", action="append", default=[], help="extra nvcc args")
    ap.add_argument("--ops", action="store_true", help="opcode histogram of the slot loops")
    a = ap.parse_args()
    defs = dict(PRESET)
    for d in a.define:
        k,_,v = d.partition("="); defs[k] = v or 1
    B = int(defs["ECC_BATCH"])
    work = tempfile.mkdtemp()
    cubin, log = build(defs, a.arch, work, a.extra)
    reg = re.search(r"Used (\d+) registers", log); spill = re.search(r"(\d+) bytes spill stores, (\d+) bytes spill loads", log)
    insts = sass(cubin, a.kernel)
    if a.dump:
        with open(a.dump,"w") as f:
            for addr, op, rest, pred in insts: f.write("%05x %s%s\n"%(addr, op, rest))
    # Subroutines: the kernel body ends at its EXIT; each callee spans CALL target .. RET.
    total = count(insts, 0, 1<<40)
    print("kernel %s: %d SASS instructions, regs %s, spills %s" % (
        a.kernel, total["n"], reg.group(1) if reg else "?", spill.groups() if spill else "?"))
    rets = [addr for addr, op, rest, pred in insts if op.startswith("RET")]
    targets = sorted(set(int(re.search(r"0x([0-9a-f]+)", rest).group(1),16)
                         for addr, op, rest, pred in insts if op.startswith("CALL")))
    subs = {}
    for t in targets:
        end = min(r for r in rets if r >= t)
        subs[t] = (t, end)
    kernel_end = min(a_ for a_, op, r_, p_ in insts if op.startswith("EXIT") and a_ > 0x1000) if targets else 1<<40
    def calls_in(lo, hi):
        c = Counter()
        for addr, op, rest, pred in insts:
            if lo <= addr <= hi and op.startswith("CALL"):
                c[int(re.search(r"0x([0-9a-f]+)", rest).group(1),16)] += 1
        return c
    # Full cost of a range = own instructions + callee cost per call site (callees may call further).
    memo = {}
    def full(lo, hi, depth=0):
        key = (lo, hi)
        if key in memo: return memo[key]
        c = count(insts, lo, hi)
        for t, n in calls_in(lo, hi).items():
            if subs[t][0] == lo: continue
            sc = full(*subs[t], depth+1)
            for k in ("n","alu","imad","clmad","mem"): c[k] += n * sc[k]
        memo[key] = c
        return c
    print("\nsubroutines (own SASS, then with nested calls):")
    print("%-14s %6s %6s %6s %6s %6s   %s" % ("range","n","alu","imad","clmad","mem","called from"))
    for t,(lo,hi) in sorted(subs.items()):
        c = count(insts, lo, hi)
        callers = [hex(addr) for addr, op, rest, pred in insts if op.startswith("CALL") and int(re.search(r"0x([0-9a-f]+)", rest).group(1),16)==t]
        print("%-6x-%-7x %6d %6d %6d %6d %6d   %s" % (lo, hi, c["n"], c["alu"], c["imad"], c["clmad"], c["mem"], " ".join(callers)))

    L = [l for l in loops(insts) if l[1] < kernel_end]
    L.sort(key=lambda x: (x[0], -x[1]))
    print("\nloops in the kernel body (own SASS):")
    print("%-10s %-10s %6s %6s %6s %6s %6s" % ("start","end","n","alu","imad","clmad","mem"))
    for lo, hi in L:
        c = count(insts, lo, hi)
        print("%-10x %-10x %6d %6d %6d %6d %6d" % (lo, hi, c["n"], c["alu"], c["imad"], c["clmad"], c["mem"]))
    step = max(L, key=lambda x: x[1]-x[0])
    inner = [l for l in L if l != step and step[0] < l[0] and l[1] < step[1]]
    inner = [l for l in inner if not any(o != l and o[0] <= l[0] and l[1] <= o[1] for o in inner)]
    inner.sort()
    cstep = full(*step)
    cin = [full(*l) for l in inner]
    rest = Counter(cstep)
    for c in cin: rest.subtract(c)
    print("\nstep loop %x-%x; slot loops: %s" % (step[0], step[1], ["%x-%x"%l for l in inner]))
    per = Counter()
    for c in cin:
        for k in ("alu","imad","clmad","mem","n"): per[k] += c[k]
    print("\nper scalar update, batch %d (callees attributed to their call sites):" % B)
    print("%-34s %8s %8s %8s %8s %8s" % ("", "n", "alu", "imad", "clmad", "mem"))
    for name, c in zip(["slot loop %d"%i for i in range(len(cin))], cin):
        print("%-34s %8d %8d %8d %8d %8d" % (name, c["n"], c["alu"], c["imad"], c["clmad"], c["mem"]))
    print("%-34s %8.1f %8.1f %8.1f %8.2f %8.1f" % ("inversion + step overhead, /B", rest["n"]/B, rest["alu"]/B, rest["imad"]/B, rest["clmad"]/B, rest["mem"]/B))
    tot = {k: per[k] + rest[k]/B for k in ("n","alu","imad","clmad","mem")}
    print("%-34s %8.1f %8.1f %8.1f %8.2f %8.1f" % ("TOTAL per update", tot["n"], tot["alu"], tot["imad"], tot["clmad"], tot["mem"]))
    print("ALU slots/update (alu + 2.01 imad) = %.1f ; clmad/update = %.2f ; instr/update = %.1f"
          % (tot["alu"]+2.01*tot["imad"], tot["clmad"], tot["n"]))
    if a.ops:
        h = Counter()
        ranges = list(inner) + [subs[t] for l in inner for t in calls_in(*l)]
        for lo, hi in ranges:
            for addr, op, rest_, pred in insts:
                if lo <= addr <= hi: h[op.split(".")[0]] += 1
        print("\nopcode histogram, slot loops + their callees (static):")
        for op, n in h.most_common(): print("  %-12s %5d" % (op, n))

if __name__ == "__main__":
    main()
