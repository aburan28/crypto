#!/usr/bin/env python3
"""Where the packed walk's static SASS cost goes, per scalar update, by source
function.  Companion to kernel_cost.py, which prices the whole kernel and its
loops; this says which routine each slot belongs to.

Builds walk() exactly as kernel_cost.py does (same preset and defines) but
with -lineinfo, disassembles with `nvdisasm -g`, and attributes:

  * every __noinline__ callee body (the products, the Frobenius networks,
    mul131) as ONE unit, named from its symbol -- the reliable form: the
    address range of a callee is exact, and it is weighted by the number of
    call sites in each slot loop, or by 1/B from the step loop, as
    kernel_cost.py weights it;
  * the kernel body's own instructions by the source function of their line
    marker, which is reliable there.

It deliberately does NOT attribute lines INSIDE callee bodies: ptxas's line
markers inside a heavily scheduled __noinline__ body smear across the header
that inlined into it (the first version of this tool put a third of
mulPolynomialPair131 under sigmaWalkNetwork131 that way).  A callee's
composition is read from its source, not from line info.

Pricing is kernel_cost.py's, plus IMNMX/VIMNMX at 1.81 slots (measured,
ITERATION-FUNCTION.md section 3.1; kernel_cost.py prices it at one) reported
as its own column and folded into the `slots` total; the kernel_cost-equivalent
total is printed for comparison.  Needs nvcc, cuobjdump and nvdisasm (CUDA
13.3+ for clmad) on PATH; see THROUGHPUT-20B.md for a toolchain that needs no
GPU and no system CUDA install.

  ./kernel_attribution.py --define ECC_WALK_TABLE=1
"""
import argparse, os, re, subprocess, sys, tempfile
from collections import Counter, defaultdict

HERE = os.environ.get("ECC2K130", os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from kernel_cost import PRESET, SPIKE, ALU, MEM, CTRL, QUARTER, QUARTER_SLOTS, IMAD_SLOTS  # noqa: E402

IMNMX_SLOTS = 1.81

def cls(op):
    h = op.split(".")[0]
    if h.startswith("CLMAD") or h.startswith("CLMUL"): return "clmad", 0.0
    if op.startswith(("IMAD.MOV", "IMAD.IADD", "IMAD.SHL")): return "alu", 1.0
    if h.startswith("IMAD") or h.startswith("IMUL"): return "imad", IMAD_SLOTS
    if h in QUARTER: return "quarter", QUARTER_SLOTS
    if h in ("IMNMX", "VIMNMX"): return "imnmx", IMNMX_SLOTS
    if h in MEM: return "mem", 0.0
    if h in CTRL: return "ctrl", 0.0
    return "alu", 1.0

def build(defs, arch, work):
    src = os.path.join(work, "walk_spike.cu"); open(src, "w").write(SPIKE)
    cubin = os.path.join(work, "walk.cubin")
    cmd = ["nvcc", "-O3", "-std=c++17", "-arch=" + arch, "-cubin", "-lineinfo", "-Xptxas", "-v",
           "-I", HERE, "-I", os.path.join(HERE, "include"), "-o", cubin, src]
    cmd[1:1] = ["-D%s=%s" % (k, v) for k, v in defs.items()]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=os.path.join(HERE, "src"))
    if r.returncode: sys.exit("nvcc failed:\n" + r.stderr[:4000])
    return cubin, r.stderr

def disassemble(cubin, kernel):
    r = subprocess.run(["nvdisasm", "-g", "-c", cubin], capture_output=True, text=True)
    if r.returncode: sys.exit(r.stderr)
    insts, cur, insec, pending, labels = [], (None, 0), False, [], {}
    for line in r.stdout.splitlines():
        if line.startswith(".text."):
            insec = kernel in line; continue
        if not insec: continue
        m = re.match(r'\s*//## File "([^"]+)", line (\d+)', line)
        if m: cur = (os.path.basename(m.group(1)), int(m.group(2))); continue
        m = re.match(r"\s*([.$A-Za-z_][\w$.]*):\s*$", line)
        if m: pending.append(m.group(1)); continue
        m = re.match(r"\s*/\*([0-9a-f]+)\*/\s+(@!?U?P\d\s+)?([A-Z][A-Z0-9_.]*)(.*)", line)
        if m:
            a = int(m.group(1), 16)
            for l in pending: labels[l] = a
            pending = []
            insts.append((a, m.group(3), m.group(4), cur))
    for l in pending: labels[l] = insts[-1][0] + 0x10
    insts = [(a, op, re.sub(r"`\(([^)]+)\)", lambda mm: "0x%x" % labels[mm.group(1)], rest), cur)
             for a, op, rest, cur in insts]
    return insts, labels

def function_map(files):
    """(file, line) -> enclosing function name, from the source text."""
    defn = re.compile(r'^\s*(?:static\s+|template\s*<[^>]*>\s*|__device__\s+|__global__\s+|__forceinline__\s+|inline\s+|ECC_HD\s+|ECC_BIG\s+|ECC_BOUNDS\s+|TW_FN\s+)*[\w:<>&*\s]+?\b([A-Za-z_]\w*)\s*\([^;{]*\)\s*(?:const\s*)?\{')
    out = {}
    for f in files:
        p = os.path.join(HERE, "include", f)
        if not os.path.exists(p): continue
        cur = "?"
        for i, line in enumerate(open(p), 1):
            m = defn.match(line)
            if m and m.group(1) not in ("if", "for", "while", "switch", "return"): cur = m.group(1)
            out[(f, i)] = cur
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arch", default="sm_120")
    ap.add_argument("--define", action="append", default=[])
    ap.add_argument("--kernel", default="4walk")
    a = ap.parse_args()
    defs = dict(PRESET)
    for d in a.define:
        k, _, v = d.partition("="); defs[k] = v or 1
    B = int(defs["ECC_BATCH"])
    work = tempfile.mkdtemp()
    cubin, log = build(defs, a.arch, work)
    reg = re.search(r"Used (\d+) registers", log)
    insts, labels = disassemble(cubin, a.kernel)
    target = lambda rest: int(re.search(r"0x([0-9a-f]+)", rest).group(1), 16)
    loops = [(target(rest), ad) for ad, op, rest, _ in insts if op.startswith("BRA") and re.search(r"0x", rest) and target(rest) <= ad]
    rets = [ad for ad, op, *_ in insts if op.startswith("RET")]
    subs = {}
    for ad, op, rest, _ in insts:
        if op.startswith("CALL"):
            t = target(rest); subs[t] = (t, min(r for r in rets if r >= t))
    kernel_end = min(ad for ad, op, *_ in insts if op.startswith("EXIT") and ad > 0x1000)
    L = sorted([l for l in loops if l[1] < kernel_end], key=lambda x: (x[0], -x[1]))
    step = max(L, key=lambda x: x[1] - x[0])
    inner = [l for l in L if l != step and step[0] < l[0] and l[1] < step[1]]
    inner = [l for l in inner if not any(o != l and o[0] <= l[0] and l[1] <= o[1] for o in inner)]
    inner.sort()
    name_at = {ad: l for l, ad in labels.items() if not l.startswith(".L_")}
    def short(sym):
        # the callee's own name is the LAST eccPacked131<len><name> in the
        # mangled symbol; the first is the kernel it was cloned into
        found = re.findall(r"eccPacked131(\d+)([A-Za-z_]\w*)", sym)
        if not found: return sym[-30:]
        n, name = found[-1]; name = name[:int(n)]
        m = re.search(r"specialized_\$_(\d+)", sym)
        return name + ("#%s" % m.group(1) if m else "")
    # weights: slot loops x1, callee bodies per call site, step remainder /B
    weight = defaultdict(float)          # kernel-body address -> weight
    callee_w = Counter()                 # callee entry -> weight
    for ad, op, rest, _ in insts:
        if ad > kernel_end: continue
        in_slot = any(lo <= ad <= hi for lo, hi in inner)
        w = 1.0 if in_slot else (1.0 / B if step[0] <= ad <= step[1] else 0.0)
        if not w: continue
        weight[ad] += w
        if op.startswith("CALL"): callee_w[target(rest)] += w
    fmap = function_map(set(c[0] for *_, c in insts if c[0]))
    rows = defaultdict(Counter)
    for ad, op, rest, (f, ln) in insts:
        if weight[ad] == 0: continue
        k, sl = cls(op); key = ("%s:%s" % (f, fmap.get((f, ln), "?")), "body")
        rows[key][k] += weight[ad]; rows[key]["slots"] += weight[ad] * sl
    for t, w in callee_w.items():
        lo, hi = subs[t]; key = (short(name_at.get(t, "?")), "callee x%.4g" % w)
        for ad, op, rest, _ in insts:
            if lo <= ad <= hi:
                k, sl = cls(op); rows[key][k] += w; rows[key]["slots"] += w * sl
    tot = Counter()
    for c in rows.values(): tot.update(c)
    print("kernel %s: regs %s; step loop %x-%x; slot loops %s; callees %d" % (
        a.kernel, reg.group(1) if reg else "?", step[0], step[1], ["%x-%x" % l for l in inner], len(subs)))
    print("%-46s %-14s %8s %8s %6s %6s %6s %7s %6s" % ("where", "kind", "slots", "alu", "imad", "quart", "imnmx", "clmad", "mem"))
    for (where, kind), c in sorted(rows.items(), key=lambda kv: -kv[1]["slots"]):
        print("%-46s %-14s %8.1f %8.1f %6.1f %6.1f %6.1f %7.2f %6.1f" % (
            where, kind, c["slots"], c["alu"], c["imad"], c["quarter"], c["imnmx"], c["clmad"], c["mem"]))
    print("%-46s %-14s %8.1f %8.1f %6.1f %6.1f %6.1f %7.2f %6.1f" % (
        "TOTAL per update", "", tot["slots"], tot["alu"], tot["imad"], tot["quarter"], tot["imnmx"], tot["clmad"], tot["mem"]))
    print("kernel_cost.py-equivalent slots (IMNMX at 1.0): %.1f" % (tot["slots"] - (IMNMX_SLOTS - 1.0) * tot["imnmx"]))

if __name__ == "__main__":
    main()
