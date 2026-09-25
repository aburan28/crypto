#!/usr/bin/env python3
"""Percent of peak for the packed walk: a per-pipe roofline from dynamic SASS counts.

The question this answers is "what fraction of the chip's throughput is the
Pollard rho kernel getting, and which unit is it short of".  kernel_cost.py and
kernel_attribution.py count static SASS weighted by loop trip counts; that
counts both arms of every branch, the distinguished-point report that almost
never runs, and a full slot-0 arm in every reverse-pass iteration.  This tool
counts what a warp actually issues per scalar update:

  1. Build the walk kernel for the requested preset.  The -D flags come from
     the Makefile itself (`make -n`), so a preset here is exactly the binary
     `make gpu-rtx-pro6000-20b` (or gpu-preset KNOBS=...) produces; the walk
     is compiled on its own, as kernel_cost.py does, which leaves its SASS
     identical to the one in the full client (checked against cuobjdump of
     the client for the 20 B/s build and with both new knobs: 4,920 and
     6,856 instructions, no difference beyond how the two disassemblers print
     branch targets).
  2. Disassemble with `nvdisasm -gi`, which gives every instruction's full
     inlining chain, and build the control-flow graph of the kernel and of
     each __noinline__ callee.
  3. Give every conditional branch a probability.  Loop back-edges get
     1 - 1/trips, the trip count read from the loop's source line (the step
     loop, the two slot loops, the fruitless-cycle retry).  Forward branches
     get the probability of the source region the branch enters: the
     distinguished-point report and the overdue restart ~0, the reverse
     pass's slot-0 arm 1/B, the pipelined forward pass's guarded tail
     (B-2)/(B-1), and so on (REGION_RULES).  Anything no rule covers is
     reported and defaults to 0.5, so an unexplained branch is visible.
  4. Solve the visit equations of each graph exactly and sum instructions x
     visits per scalar update.  Predicated-off instructions still issue, so
     they count in full.  A self-check: the dynamic CLMAD count must come out
     at the value the arithmetic implies (33.125 per update for the 20 B/s
     build: 77 products + 8 inversion products per 16 updates, 6 CLMADs
     each, plus 20 inversion squarings).
  5. Price every instruction on the pipe that executes it, at the measured
     lane rate of that pipe on the target GPU (MACHINES, with the receipt
     each number comes from), and report SM-clocks per update per pipe, the
     ceiling that pipe alone imposes, and -- given a measured rate -- how
     busy each pipe is.

Ceilings are upper bounds from counting, not predictions: a kernel can sit
below every ceiling because of latency and dependencies.  The utilisation
column is the useful output: a pipe at ~100% is the binding one, and the
largest utilisation says how much any other change can buy before that pipe
is the wall.  Where a pipe's peak is not known to one number the machine
model carries the range the receipts support (the carry-less unit on the RTX
PRO 6000 has been measured at 1.62, 1.69 and 1.99 lanes per SM-clock
depending on the instruction pattern), and the report prints every one.

    ./roofline.py                                  # 20 B/s build, RTX PRO 6000
    ./roofline.py --measured 20.078                # with the receipt's rate
    ./roofline.py --knobs "PACKED_ALU_SQR=1"       # any gpu-preset override
    ./roofline.py --target gpu-b200-19b --gpu b200 --measured 19.40
    ./roofline.py --make-var PACKED_SQUARE_TABLE=1 --make-var PACKED_INV_POLY=1
    ./roofline_calibrate.py                        # every frozen sweep receipt against each pipe model
    ./roofline.py --json out.json                  # machine-readable

Needs nvcc, nvdisasm (CUDA 13.3+ for clmad) on PATH; no GPU.
"""
import argparse
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
KERNEL = "_ZN12eccPacked1314walkE10WalkParamsIjEPj"
SPIKE = '#include "../include/curveparams.h"\n#include "../include/packedkernels.cuh"\n'

# --------------------------------------------------------------------------
# Machine models.  Rates are lanes per SM-clock; every number names the
# receipt it comes from.  "issue" is 4 schedulers x 32 lanes.
# --------------------------------------------------------------------------
MACHINES = {
    "rtx-pro-6000": {
        "name": "RTX PRO 6000 Blackwell Server Edition (sm_120)",
        "arch": "sm_120",
        "sms": 188,
        # SM clock sampled after every timed run of the 20 B/s build in
        # benchmarks/two-chains/summary.json and autosweep/rtx-pro-6000: 2400-2422.
        "clock_ghz": 2.415,
        "clock_source": "benchmarks/two-chains/summary.json ref smClockMHz 2400-2422 (median 2415)",
        "pipes": {
            "issue": {"rate": 128.0, "source": "4 schedulers x 1 warp-instruction per clock"},
            "alu": {"rate": 64.0, "source": "benchmarks/hardware-limits/result.json: lop3 62.2, iadd3 63.0, "
                    "shf 62.5, prmt 62.2 lanes/SM-clk; benchmarks/fast-clmad/probe LOP3 63.9"},
            "fma": {"rate": 64.0, "source": "benchmarks/hardware-limits/result.json: imad 62.3 lanes/SM-clk"},
            "xu": {"rate": 16.0, "source": "benchmarks/fast-clmad/probe POPC 16.1, FLO 15.9 lanes/SM-clk"},
            "lsu": {"rate": 16.0, "source": "benchmarks/hardware-limits/result.json ld.shared 15.9 lane-loads/SM-clk"},
            "clmad": {"rate": 1.99, "source": "benchmarks/fast-clmad/probe CLMAD lo+hi product stream 1.99",
                      "alternatives": {
                          "1.62 (ONE-BLOCK-GEOMETRY.md section 1, 'kernel mix')": 1.62,
                          "1.69 (fast-clmad probe, CLMAD.lo-only stream)": 1.69,
                          "1.99 (fast-clmad probe, lo+hi product stream)": 1.99,
                          "2.00 (FP64 pipe at 1/64 of FP32: 2 lanes/SM-clk)": 2.00}},
        },
    },
    "b200": {
        "name": "B200 (sm_100)",
        "arch": "sm_100",
        "sms": 148,
        "clock_ghz": 1.965,
        "clock_source": "benchmarks/autosweep/b200/summary.json smClockMHz 1965 throughout",
        "pipes": {
            "issue": {"rate": 128.0, "source": "4 schedulers x 1 warp-instruction per clock"},
            "alu": {"rate": 63.7, "source": "benchmarks/fast-clmad/probe LOP3 63.7 lanes/SM-clk"},
            "fma": {"rate": 64.0, "source": "assumed as sm_120 (IMAD.WIDE 22.9 in the probe is the wide form)"},
            "xu": {"rate": 16.0, "source": "benchmarks/fast-clmad/probe POPC 16.0"},
            "lsu": {"rate": 16.0, "source": "assumed as sm_120"},
            "clmad": {"rate": 29.1, "source": "benchmarks/fast-clmad/probe CLMAD lo+hi product stream 29.1",
                      "alternatives": {"16.9 (CLMAD.lo-only stream)": 16.9, "29.1 (lo+hi product stream)": 29.1}},
        },
    },
}

# SASS opcode -> (pipe, lane-slots per lane-instruction).  The slot weight is the
# pipe's single-op rate over this op's measured rate (IMNMX 35.4 against LOP3's
# 63.9 in the fast-clmad probe; the .WIDE/.HI forms of IMAD are two passes).
def pipe_of(op):
    head = op.split(".")[0]
    mods = op.split(".")[1:]
    if head in ("CLMAD", "CLMUL"):
        return "clmad", 1.0
    if head in ("IMAD", "IMUL", "IMADSP", "IMNMX32I") or head.startswith("IMAD"):
        if "WIDE" in mods or "HI" in mods:
            return "fma", 2.0
        return "fma", 1.0
    if head in ("FFMA", "FADD", "FMUL", "FFMA32I", "FADD32I", "FMUL32I", "HFMA2", "HADD2", "HMUL2",
                "FSWZADD", "FCHK"):
        return "fma", 1.0
    if head in ("DFMA", "DADD", "DMUL", "DSETP"):
        return "clmad", 1.0  # the FP64 pipe that carries CLMAD
    if head in ("POPC", "FLO", "BREV", "MUFU", "I2F", "F2I", "I2FP", "F2F", "FRND", "I2I", "F2IP"):
        return "xu", 1.0
    if head in ("IMNMX", "VIMNMX", "VIMNMX3"):
        return "alu", 1.8
    if head in ("LDS", "STS", "LDG", "STG", "LDL", "STL", "LD", "ST", "ATOM", "ATOMS", "ATOMG", "RED",
                "REDG", "SHFL", "LDSM", "LDGSTS", "MEMBAR", "CCTL", "ERRBAR", "LDGDEPBAR", "MATCH", "QSPC"):
        return "lsu", 1.0
    if head in ("LDC", "LDCU"):
        return "const", 1.0
    if head.startswith("U") or head in ("S2UR", "R2UR", "VOTEU", "CS2R", "S2R"):
        # uniform datapath and special-register reads: an issue slot, no vector pipe
        return "uniform", 1.0
    if head in ("BRA", "BRX", "JMP", "JMX", "CALL", "RET", "EXIT", "BSSY", "BSYNC", "WARPSYNC", "BAR",
                "BMOV", "BPT", "YIELD", "NANOSLEEP", "DEPBAR", "KILL", "BREAK", "SYNCS", "ACQBULK"):
        return "control", 1.0
    if head == "NOP":
        return "nop", 1.0
    # LOP3, SHF, PRMT, IADD3, ISETP, LEA, SEL, MOV, BMSK, SGXT, PLOP3, P2R, R2P, VOTE, ...
    return "alu", 1.0


# --------------------------------------------------------------------------
# Build
# --------------------------------------------------------------------------
def tool(name):
    return os.environ.get(name.upper(), name)


def make_defines(target, knobs, arch_override=None, make_vars=None):
    """The exact -D flags `make <target> KNOBS=...` compiles the client with."""
    cmd = ["make", "-n", "-s", target] + list(make_vars or [])
    if knobs:
        cmd.append("KNOBS=" + knobs)
    if arch_override:
        cmd.append("PRO6000_ARCH=" + arch_override)
    first = subprocess.run(cmd, cwd=HERE, capture_output=True, text=True)
    if first.returncode:
        sys.exit("make -n %s failed:\n%s" % (target, first.stderr))
    joined = re.sub(r"\\\n\s*", " ", first.stdout)   # recipe lines continued with a backslash
    line = next((l for l in joined.splitlines() if re.search(r"\bmake\b.*\becc2k130\b", l)), None)
    if line is None:
        sys.exit("could not find the recursive make line for %s:\n%s" % (target, first.stdout))
    args = shlex.split(line)
    args = args[args.index("ecc2k130") + 1:] if "ecc2k130" in args else args
    # The child make of a real build also inherits the parent's command-line
    # variables through MAKEFLAGS, below the ones its recipe spells out.
    args = list(make_vars or []) + args
    second = subprocess.run(["make", "-n", "-B", "ecc2k130"] + args, cwd=HERE, capture_output=True, text=True)
    if second.returncode:
        sys.exit("make -n ecc2k130 failed:\n" + second.stderr)
    nv = " ".join(l.strip().rstrip("\\") for l in second.stdout.splitlines() if "nvcc" in l or "-DECC" in l)
    defs = re.findall(r"-D(ECC_[A-Z0-9_]+)=(\S+)", nv)
    arch = re.search(r"code=(sm_\d+a?)", nv)
    return dict(defs), (arch.group(1) if arch else "sm_120")


def build(defs, arch, work):
    src = os.path.join(work, "walk_spike.cu")
    with open(src, "w") as f:
        f.write(SPIKE)
    cubin = os.path.join(work, "walk.cubin")
    cmd = [tool("nvcc"), "-O3", "-std=c++17", "-arch=" + arch, "-cubin", "-lineinfo", "-Xptxas", "-v",
           "-I", HERE, "-I", os.path.join(HERE, "include"), "-o", cubin, src]
    cmd[1:1] = ["-D%s=%s" % kv for kv in sorted(defs.items())]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=os.path.join(HERE, "src"))
    if r.returncode:
        sys.exit("nvcc failed:\n" + r.stderr[:4000])
    log = r.stderr
    m = re.search(r"Compiling entry function '%s'.*?Used (\d+) registers" % re.escape(KERNEL), log, re.S)
    regs = int(m.group(1)) if m else None
    m = re.search(r"Compiling entry function '%s'.*?(\d+) bytes spill stores, (\d+) bytes spill loads"
                  % re.escape(KERNEL), log, re.S)
    spills = (int(m.group(1)), int(m.group(2))) if m else None
    dis = subprocess.run([tool("nvdisasm"), "-gi", "-c", cubin], capture_output=True, text=True)
    if dis.returncode:
        sys.exit("nvdisasm failed:\n" + dis.stderr[:2000])
    return dis.stdout, regs, spills


# --------------------------------------------------------------------------
# Parse nvdisasm -gi output
# --------------------------------------------------------------------------
class Insn:
    __slots__ = ("addr", "pred", "op", "args", "chain", "labels", "target", "fn")

    def __init__(self, addr, pred, op, args, chain, labels, fn):
        self.addr, self.pred, self.op, self.args = addr, pred, op, args
        self.chain, self.labels, self.fn = chain, labels, fn
        m = re.search(r"`\(([^)]+)\)", args)
        self.target = m.group(1) if m else None


FRAME = re.compile(r'File "([^"]+)", line (\d+)')
INSN = re.compile(r"^\s*/\*([0-9a-f]+)\*/\s+(@!?U?P[T0-9]\s+)?([A-Z][A-Z0-9_.]*)\s*(.*?)\s*;?\s*$")


def parse(text, section):
    """Functions of one .text section, in address order: {name: [Insn]}."""
    funcs, order = {}, []
    cur_fn, in_section, chain, labels, last_chain = None, False, [], [], ()
    for line in text.splitlines():
        s = line.strip()
        if s.startswith(".text."):
            in_section = s[len(".text."):].rstrip(":") == section
            cur_fn = None
            continue
        if not in_section:
            continue
        if s.startswith("//##"):
            m = FRAME.search(s)
            if m:
                chain.append((os.path.realpath(m.group(1)), int(m.group(2))))
            continue
        if s.endswith(":") and not s.startswith("/*") and " " not in s:
            name = s[:-1]
            if name.startswith(".L_"):
                labels.append(name)
            else:
                cur_fn = name
                funcs[cur_fn] = []
                order.append(cur_fn)
                labels, last_chain = [name], ()
            continue
        m = INSN.match(line)
        if m and cur_fn is not None:
            pred = (m.group(2) or "").strip()
            if chain:
                last_chain = tuple(chain)
            funcs[cur_fn].append(Insn(int(m.group(1), 16), pred, m.group(3), m.group(4), last_chain,
                                      tuple(labels), cur_fn))
            chain, labels = [], []
    return funcs, order


# --------------------------------------------------------------------------
# Source regions: where a branch goes, and how often
# --------------------------------------------------------------------------
_src_cache = {}


def src_lines(path):
    if path not in _src_cache:
        try:
            with open(path) as f:
                _src_cache[path] = f.read().split("\n")
        except OSError:
            _src_cache[path] = []
    return _src_cache[path]


def block_extent(lines, i, col):
    """1-based [first, last] of the construct whose condition starts at lines[i][col:].
    A braced body runs to its matching brace; an unbraced one to its ';'."""
    depth, j, c, seen_paren = 0, i, col, False
    # skip the parenthesised condition, if any
    while j < len(lines):
        L = lines[j]
        while c < len(L):
            ch = L[c]
            if ch == "(":
                depth += 1
                seen_paren = True
            elif ch == ")":
                depth -= 1
                if depth == 0 and seen_paren:
                    c += 1
                    break
            c += 1
        else:
            j, c = j + 1, 0
            continue
        break
    # body: '{' ... matching '}', or a statement to ';'
    while j < len(lines):
        L = lines[j]
        while c < len(L) and L[c] in " \t":
            c += 1
        if c < len(L):
            break
        j, c = j + 1, 0
    if j >= len(lines):
        return None
    if lines[j][c] == "{":
        depth = 0
        while j < len(lines):
            L = lines[j]
            while c < len(L):
                if L[c] == "{":
                    depth += 1
                elif L[c] == "}":
                    depth -= 1
                    if depth == 0:
                        return (i + 1, j + 1)
                c += 1
            j, c = j + 1, 0
        return None
    depth = 0
    while j < len(lines):
        L = lines[j]
        while c < len(L):
            if L[c] in "({":
                depth += 1
            elif L[c] in ")}":
                depth -= 1
            elif L[c] == ";" and depth == 0:
                return (i + 1, j + 1)
            c += 1
        j, c = j + 1, 0
    return None


def body_extent(lines, ext):
    """Narrow a braced construct's extent to its body: the lines after the one
    holding the opening brace, through the closing brace."""
    lo, hi = ext
    # A closing line that goes on ("} else if (...) {", "} else {") holds the next
    # construct's condition, not this body's code.
    if lines[hi - 1].strip().startswith("}") and lines[hi - 1].strip() != "}" and hi > lo:
        hi -= 1
    for j in range(lo - 1, hi):
        L = lines[j]
        k = L.find("{")
        if k >= 0:
            return (j + 2, hi) if not L[k + 1:].strip() else (j + 1, hi)
    return (lo, hi)


def find_regions(rules, batch, steps=1024):
    """[(path, lo, hi, prob, name)] for every occurrence of every rule's pattern."""
    out = []
    for fname, pattern, prob_expr, name in rules:
        path = os.path.realpath(os.path.join(HERE, "include", fname))
        lines = src_lines(path)
        for i, L in enumerate(lines):
            k = L.find(pattern)
            if k < 0:
                continue
            ext = block_extent(lines, i, k)
            if ext is None:
                continue
            lo, hi = body_extent(lines, ext)
            prob = eval(prob_expr, {"B": batch, "U": 1, "S": steps, "P_RETRY": P_RETRY})
            out.append((path, lo, hi, prob, name, prob_expr))
    # An else body is entered when its if is not: find the if each "} else {"
    # closes and give its body the complement of that if's probability.
    for path in sorted({r[0] for r in out}):
        lines = src_lines(path)
        for i, L in enumerate(lines):
            if not L.strip().startswith("} else {"):
                continue
            depth, j = 0, i
            col = L.index("}")
            opener = None
            while j >= 0 and opener is None:
                text = lines[j] if j != i else lines[j][:col + 1]
                for c in reversed(text):
                    if c == "}":
                        depth += 1
                    elif c == "{":
                        depth -= 1
                        if depth == 0:
                            opener = j
                            break
                j -= 1
            if opener is None:
                continue
            rule = next((r for r in out if r[0] == path and r[1] <= opener + 2 <= r[2] + 1
                         and r[1] - 1 <= opener + 1), None)
            rule = next((r for r in out if r[0] == path and r[1] in (opener + 1, opener + 2)
                         and r[4] != "live walk"), rule)
            if rule is None:
                continue
            ext = block_extent(lines, i, L.index("else") + 4)
            if ext is None:
                continue
            lo, hi = body_extent(lines, (i + 1, ext[1]))
            expr = "1.0-(%s)" % rule[5]
            prob = eval(expr, {"B": batch, "U": 1, "S": steps, "P_RETRY": P_RETRY})
            out.append((path, i + 2, hi, prob, "else of: " + rule[4], expr))
    return out


# (file, pattern, probability of entering the construct from its condition, name).
# B is ECC_BATCH.  Probabilities are per warp: a region any lane enters costs the
# whole warp, so "~0" means rare for every one of 32 lanes.
REGION_RULES = [
    ("packedkernels.cuh", "if (hw <= p.dpWeight)", "0.0", "distinguished-point report (~2^-23 per warp-slot)"),
    ("packedkernels.cuh", "} else if (guard && now - p.startIter[id] >= p.maxIters)", "0.0",
     "overdue restart (once per ECC_GUARD_PERIOD steps, then rare)"),
    ("packedkernels.cuh", "if (!p.dead[id])", "1.0", "live walk"),
    ("packedkernels.cuh", "if (slot > 1) store(p.pchain, slot - 1", "(B-2)/(B-1)", "pipelined W store"),
    ("packedkernels.cuh", "if (slot + 1 < ECC_BATCH)", "(B-2)/(B-1)", "pipelined next selection"),
    # U is the unroll of the enclosing slot loop: only the body copy that meets
    # slot 0 keeps the branch, and it meets it once in B/U iterations.
    ("packedkernels.cuh", "if (slot) {", "(B-U)/B", "slot > 0 arm"),
    ("packedkernels.cuh", "if (slot > 0) {", "(B-U)/B", "slot > 0 prefetch"),
    # the two-chain kernel (PACKED_CHAINS=2): L = B/2 slots per chain
    ("packedkernels.cuh", "if (i) {", "(B/2.0-U)/(B/2.0)", "chains: slot > 0 arm"),
    ("packedkernels.cuh", "if (i + 1 < L)", "(B/2.0-2)/(B/2.0-1)", "chains: next selection"),
    ("packedkernels.cuh", "if (!first) {", "(B-1)/B", "fused: not the first slot"),
    ("packedkernels.cuh", "if (i + 1 < ECC_BATCH) {", "(B-1)/B", "fused: not the last slot"),
    ("packedkernels.cuh", "if (!last)", "(S-1.0)/S", "fused: not the launch's last step"),
    ("packedkernels.cuh", "if (!slot) {", "1.0*U/B", "slot 0 pipeline refresh"),
    ("packedtablewalk.cuh", "for (int i = 0; i < TW_H && eccTagFruitless(tag, win, rpow, 131, TW_RM); ++i)", "P_RETRY", "fruitless-cycle retry"),
]
# A step is refused when it undoes the last one (1 in 2*H*m) or closes a
# longer fruitless run (rarer; spurious refusals, 3 in 2^16, rarer still);
# WALK-CONSTANT.md measures the whole rule.  The retry loop
# is priced at 1/(2*8*131); any lane retrying costs the warp one more pass.
P_RETRY = 1.0 / (2 * 8 * 131)

# Loops: (file, pattern of the loop statement, trip count per entry).  A back-edge
# belongs to the innermost of these constructs that holds most of the
# instructions between its target and itself.
LOOP_RULES = [
    ("packedkernels.cuh", "for (int step = 0; step < p.steps; ++step)", "STEPS"),
    ("packedkernels.cuh", "for (int slot = 1; slot < ECC_BATCH; ++slot)", "B-1"),
    ("packedkernels.cuh", "for (int slot = 0; slot < ECC_BATCH; ++slot)", "B"),
    ("packedkernels.cuh", "for (int slot = ECC_BATCH - 1; slot >= 0; --slot)", "B"),
    ("packedkernels.cuh", "for (int i = 0; i < ECC_BATCH; ++i)", "B"),
    ("packedkernels.cuh", "for (int i = 1; i < L; ++i)", "B/2-1"),
    ("packedkernels.cuh", "for (int i = L - 1; i >= 0; --i)", "B/2"),
    ("packedkernels.cuh", "for (int slot = ECC_BATCH / 2 - 1; slot >= 0; --slot)", "B/2"),
    ("packedkernels.cuh", "for (int k = 1; k < ECC_BATCH / 2; ++k)", "B/2-1"),
    ("packedtablewalk.cuh", "for (int i = 0; i < TW_H && eccTagFruitless(tag, win, rpow, 131, TW_RM); ++i)", "RETRY"),
    ("packedtablewalk.cuh", "for (int i = threadIdx.x; i < words; i += blockDim.x)", "SMEM_FILL"),
]


def unroll_factor(lines, i, unroll_slots):
    """The unroll the pragma block above loop line i asks for: an explicit
    `#pragma unroll N`, or ECC_UNROLL_SLOTS where the block selects on it."""
    for j in range(i - 1, max(-1, i - 12), -1):
        t = lines[j].strip()
        if not t.startswith("#"):
            break
        if "ECC_UNROLL_SLOTS" in t:
            return max(1, min(8, 1 << (max(1, unroll_slots).bit_length() - 1)))
        m = re.match(r"#pragma unroll (\d+)", t)
        if m and j == i - 1:
            return int(m.group(1))
    return 1


def find_loops(batch, steps, threads, tw_words, unroll_slots=1):
    out = []
    for fname, pattern, expr in LOOP_RULES:
        path = os.path.realpath(os.path.join(HERE, "include", fname))
        lines = src_lines(path)
        for i, L in enumerate(lines):
            k = L.find(pattern)
            if k < 0:
                continue
            ext = block_extent(lines, i, k)
            if ext is None:
                continue
            B = batch
            trips = {"STEPS": steps, "B-1": B - 1, "B": B, "B/2": B // 2, "B/2-1": B // 2 - 1,
                     "RETRY": None, "SMEM_FILL": max(1, -(-tw_words // threads))}[expr]
            u = unroll_factor(lines, i, unroll_slots) if expr.startswith("B") else 1
            if trips and u > 1:
                trips = max(1, trips // u)
            out.append((path, ext[0], ext[1], trips, pattern + ("  (unrolled x%d)" % u if u > 1 else "")))
    return out


def frame_text(frame):
    path, line = frame
    lines = src_lines(path)
    return lines[line - 1] if 0 < line <= len(lines) else ""


# --------------------------------------------------------------------------
# Control-flow graph and exact visit counts
# --------------------------------------------------------------------------
def blocks_of(insns):
    """Split a function into basic blocks: [(start, end)] indices, end exclusive."""
    leaders = {0}
    for i, x in enumerate(insns):
        if x.labels and i:
            leaders.add(i)
        head = x.op.split(".")[0]
        if head in ("BRA", "BRX", "JMP", "JMX", "RET", "EXIT") and i + 1 < len(insns):
            leaders.add(i + 1)
    starts = sorted(leaders)
    return [(s, starts[k + 1] if k + 1 < len(starts) else len(insns)) for k, s in enumerate(starts)]


def is_conditional(x):
    if x.pred and x.pred not in ("@PT",):
        return True
    # BRA.U !UP0, `(target) : a uniform-predicate branch
    return bool(re.match(r"^!?U?P[0-9T]\s*,", x.args))


def solve(n, edges, entry):
    """Expected visits v = e_entry + P^T v, by sparse Gaussian elimination in floats."""
    rows = [defaultdict(float) for _ in range(n)]
    for i in range(n):
        rows[i][i] += 1.0
    for (u, w), p in edges.items():
        rows[w][u] -= float(p)
    b = [0.0] * n
    b[entry] = 1.0
    colrows = defaultdict(set)
    for r in range(n):
        for c in rows[r]:
            colrows[c].add(r)
    for col in range(n):
        piv = rows[col].get(col, 0.0)
        if abs(piv) < 1e-300:
            continue
        prow = rows[col]
        for r in list(colrows[col]):
            if r == col:
                continue
            f = rows[r].get(col, 0.0) / piv
            if f == 0.0:
                continue
            rr = rows[r]
            for c, val in prow.items():
                nv = rr.get(c, 0.0) - f * val
                if abs(nv) < 1e-15:
                    if c in rr:
                        del rr[c]
                        colrows[c].discard(r)
                else:
                    if c not in rr:
                        colrows[c].add(r)
                    rr[c] = nv
            b[r] -= f * b[col]
    return [b[i] / rows[i][i] if rows[i].get(i, 0.0) else 0.0 for i in range(n)]


class Model:
    def __init__(self, funcs, batch, steps, threads, tw_words, unroll_slots=1):
        self.funcs, self.B, self.steps = funcs, batch, steps
        self.regions = find_regions(REGION_RULES, batch, steps)
        self.threads, self.tw_words = threads, tw_words
        self.loop_regions = find_loops(batch, steps, threads, tw_words, unroll_slots)
        self.unknown, self.loops, self.branches = [], [], []

    def region_set(self, chain):
        out = set()
        for path, line in chain:
            for k, (rp, lo, hi, prob, name, expr) in enumerate(self.regions):
                if rp == path and lo <= line <= hi:
                    out.add(k)
        return out

    def block_regions(self, insns, s, e):
        """Regions of a block: the most common region set among its instructions,
        counting only chains that reach packedkernels.cuh when any do (ptxas
        truncates the inline chain of merged instructions to a lone frame)."""
        body = [x for x in insns[s:e] if x.chain]
        full = [x for x in body if any(p.endswith("packedkernels.cuh") for p, _ in x.chain)]
        c = Counter(frozenset(self.region_set(x.chain)) for x in (full or body))
        return set(c.most_common(1)[0][0]) if c else set()

    GLUE = ("BRA", "BSYNC", "BSSY", "NOP", "WARPSYNC", "JMP")

    def succ_regions(self, insns, blocks, label_at, k, depth=0):
        """Regions of the code a branch successor leads into.  A glue block --
        only jumps and reconvergence -- carries whatever source line ptxas left
        on it, so look through it: to its target if it jumps unconditionally,
        to what both sides share if it branches."""
        s, e = blocks[k]
        if depth < 6 and all(x.op.split(".")[0] in self.GLUE for x in insns[s:e]):
            last = insns[e - 1]
            nxt = k + 1 if k + 1 < len(blocks) else None
            head = last.op.split(".")[0]
            if head in ("BRA", "JMP") and last.target in label_at:
                t = label_at[last.target]
                if not is_conditional(last):
                    return self.succ_regions(insns, blocks, label_at, t, depth + 1)
                a = self.succ_regions(insns, blocks, label_at, t, depth + 1)
                b = self.succ_regions(insns, blocks, label_at, nxt, depth + 1) if nxt is not None else a
                return a & b
            if nxt is not None:
                return self.succ_regions(insns, blocks, label_at, nxt, depth + 1)
        return self.block_regions(insns, s, e)

    def region_prob(self, ks, unroll=1):
        p = Fraction(1)
        for k in ks:
            expr = self.regions[k][5]
            val = eval(expr, {"B": self.B, "U": unroll, "S": self.steps, "P_RETRY": P_RETRY})
            p *= Fraction(val).limit_denominator(10**7)
        return p

    def unroll_at(self, chain):
        """Unroll factor of the innermost loop construct holding this chain."""
        best = None
        for path, line in chain:
            for lp, lo, hi, n, pat in self.loop_regions:
                if lp == path and lo <= line <= hi and (best is None or hi - lo < best[0]):
                    m = re.search(r"unrolled x(\d+)", pat)
                    best = (hi - lo, int(m.group(1)) if m else 1)
        return best[1] if best else 1

    def loop_hits(self, chain):
        hit = set()
        for path, line in chain:
            for k, (lp, lo, hi, n, pat) in enumerate(self.loop_regions):
                if lp == path and lo <= line <= hi:
                    hit.add(k)
        return hit

    def assign_loop(self, body, backedge, claimed):
        """The source loop a SASS loop implements.  Candidates are the loop
        constructs holding at least 90% of the spanned instructions that lie in
        any loop construct (an outer loop's body is mostly its inner loops, so a
        majority rule would pick those), less the ones enclosing SASS loops have
        already claimed; the innermost wins.  With none left -- a loop inlined
        into another, whose instructions all sit in the outer construct too --
        the innermost unclaimed construct the back-edge itself lies in.  Chains
        ptxas truncated to a lone helper frame do not vote."""
        votes, voters = Counter(), 0
        for x in body:
            hit = self.loop_hits(x.chain)
            if hit:
                voters += 1
            for k in hit:
                votes[k] += 1
        span = lambda k: self.loop_regions[k][2] - self.loop_regions[k][1]
        cands = [k for k, c in votes.items() if 10 * c >= 9 * voters and k not in claimed]
        if not cands:
            cands = [k for k in self.loop_hits(backedge.chain) if k not in claimed]
        if not cands:
            cands = [k for k in votes if k not in claimed]
        return min(cands, key=span) if cands else None

    def visits(self, fname):
        insns = self.funcs[fname]
        blocks = blocks_of(insns)
        index = {s: k for k, (s, e) in enumerate(blocks)}
        label_at = {}
        for k, (s, e) in enumerate(blocks):
            for lab in insns[s].labels:
                label_at[lab] = k
        # SASS loops (back-edges), matched to source loops outermost first so an
        # inner loop cannot take a construct its enclosing loop implements.
        back = []
        for k, (s, e) in enumerate(blocks):
            last = insns[e - 1]
            if last.op.split(".")[0] in ("BRA", "JMP") and last.target in label_at:
                t = label_at[last.target]
                if insns[blocks[t][0]].addr <= last.addr:
                    back.append((k, t, blocks[t][0], e))
        assigned = {}
        for k, t, s0, e0 in sorted(back, key=lambda b: -(b[3] - b[2])):
            claimed = {assigned[k2] for k2, t2, s2, e2 in back
                       if k2 in assigned and s2 <= s0 and e0 <= e2 and (s2, e2) != (s0, e0)}
            assigned[k] = self.assign_loop(insns[s0:e0], insns[e0 - 1], claimed)
        # A loop closed by an unconditional back-edge leaves through a forward
        # branch out of its range: that branch is taken once per 1/trips.
        exit_prob = {}
        for k, t, s0, e0 in back:
            if is_conditional(insns[e0 - 1]) or assigned.get(k) is None:
                continue
            n = self.loop_regions[assigned[k]][3]
            lo, hi = insns[s0].addr, insns[e0 - 1].addr
            exits = []
            for k2, (s2, e2) in enumerate(blocks):
                x = insns[e2 - 1]
                if not (lo <= x.addr <= hi) or k2 == k:
                    continue
                if x.op.split(".")[0] in ("BRA", "JMP") and is_conditional(x) and x.target in label_at:
                    ta = insns[blocks[label_at[x.target]][0]].addr
                    if ta > hi or ta < lo:
                        exits.append(k2)
                elif x.op.split(".")[0] == "EXIT" and is_conditional(x):
                    exits.append(k2)        # a loop that returns when it is done
            if exits:
                # the loop test is the last exit before the back-edge
                # a retry loop (no trip count) leaves unless the retry fires
                exit_prob[max(exits)] = (1 - Fraction(P_RETRY).limit_denominator(10**6) if n is None
                                         else Fraction(1, max(1, n)))
                self.loops.append((fname, lo, hi, self.loop_regions[assigned[k]][4] + "  (unconditional back-edge)", n))
        edges = defaultdict(Fraction)
        for k, (s, e) in enumerate(blocks):
            last = insns[e - 1]
            head = last.op.split(".")[0]
            nxt = k + 1 if k + 1 < len(blocks) else None
            if head in ("RET",) or (head == "EXIT" and not is_conditional(last)):
                continue
            if head == "EXIT":  # predicated exit: a loop's end, else inactive threads only
                if nxt is not None:
                    edges[(k, nxt)] += 1 - exit_prob.get(k, 0)
                continue
            if head in ("BRA", "JMP") and last.target in label_at:
                t = label_at[last.target]
                if not is_conditional(last):
                    edges[(k, t)] += 1
                    continue
                # conditional
                if insns[blocks[t][0]].addr > last.addr and k in exit_prob:
                    p = exit_prob[k]
                    edges[(k, t)] += p
                    if nxt is not None:
                        edges[(k, nxt)] += 1 - p
                    continue
                if insns[blocks[t][0]].addr <= last.addr:
                    li = assigned.get(k)
                    pat, n = (self.loop_regions[li][4], self.loop_regions[li][3]) if li is not None else (None, None)
                    if pat is None and not any(x.chain for x in insns[blocks[t][0]:e]):
                        pat, n = "prologue loop without line info", 1
                    if pat is None:
                        p = Fraction(9, 10)
                        self.unknown.append((fname, last.addr, "back-edge, no loop rule", "0.9"))
                    elif n is None:
                        p = Fraction(P_RETRY).limit_denominator(10**6)
                    else:
                        p = Fraction(n - 1, n) if n > 0 else Fraction(0)
                    self.loops.append((fname, insns[blocks[t][0]].addr, last.addr, pat, n))
                    edges[(k, t)] += p
                    if nxt is not None:
                        edges[(k, nxt)] += 1 - p
                    continue
                # forward conditional: which side enters a region the other does not?
                RT = self.succ_regions(insns, blocks, label_at, t)
                RF = self.succ_regions(insns, blocks, label_at, nxt) if nxt is not None else set()
                rt, rf = RT - RF, RF - RT
                u = self.unroll_at(last.chain)
                pt, pf = self.region_prob(rt, u) if rt else None, self.region_prob(rf, u) if rf else None
                if pt is not None and pf is None:
                    p = pt
                elif pf is not None and pt is None:
                    p = 1 - pf
                elif pt is not None and pf is not None:
                    p = pt / (pt + pf) if pt + pf else Fraction(1, 2)
                else:
                    p = Fraction(1, 2)
                    self.unknown.append((fname, last.addr, "forward branch enters no known region",
                                         frame_text(last.chain[0]).strip()[:70] if last.chain else ""))
                self.branches.append((fname, last.addr, float(p),
                                      [self.regions[i][4] for i in rt], [self.regions[i][4] for i in rf]))
                edges[(k, t)] += p
                if nxt is not None:
                    edges[(k, nxt)] += 1 - p
                continue
            if nxt is not None:
                edges[(k, nxt)] += 1
        v = solve(len(blocks), edges, 0)
        per_insn = []
        for k, (s, e) in enumerate(blocks):
            for x in insns[s:e]:
                per_insn.append((x, v[k]))
        return per_insn


def dynamic_counts(funcs, kernel, batch, steps, threads, tw_words, unroll_slots=1):
    """[(Insn, executions per scalar update)] over the kernel and every callee it reaches."""
    m = Model(funcs, batch, steps, threads, tw_words, unroll_slots)
    per_call = {}

    def expand(fname, scale, out, depth=0):
        if depth > 8:
            return
        if fname not in per_call:
            per_call[fname] = m.visits(fname)
        for x, v in per_call[fname]:
            w = float(v) * scale
            out.append((x, w))
            if x.op.startswith("CALL") and x.target:
                callee = x.target
                if callee in funcs:
                    expand(callee, w, out, depth + 1)

    out = []
    expand(kernel, 1.0 / (steps * batch), out)
    return out, m


# --------------------------------------------------------------------------
# Attribution to source functions (innermost frame)
# --------------------------------------------------------------------------
FUNC_DEF = re.compile(r"^\s*(?:template\s*<[^>]*>\s*)?(?:static\s+|inline\s+|__device__\s+|__forceinline__\s+|"
                      r"__global__\s+|ECC_HD\s+|ECC_BIG\s+|TW_FN\s+|ECC_POLY_SINGLE\s+|ECC_POLY_PAIR\s+|"
                      r"const\s+|unsigned\s+|signed\s+)*[A-Za-z_][\w:<>]*[\s*&]+([A-Za-z_]\w*)\s*\([^;]*$")
FUNC_DEF1 = re.compile(r"^\s*(?:template\s*<[^>]*>\s*)?(?:(?:static|inline|__device__|__forceinline__|__global__|"
                       r"ECC_HD|ECC_BIG|TW_FN|const|unsigned|signed)\s+)*[A-Za-z_][\w:<>]*[\s*&]+"
                       r"([A-Za-z_]\w*)\s*\([^;{]*\)\s*(?:const\s*)?\{")
_fn_cache = {}


def function_at(path, line):
    if path not in _fn_cache:
        names = []
        for i, L in enumerate(src_lines(path)):
            m = FUNC_DEF.match(L) or FUNC_DEF1.match(L)
            if m and m.group(1) not in ("if", "for", "while", "switch", "return", "sizeof"):
                names.append((i + 1, m.group(1)))
        _fn_cache[path] = names
    best = "?"
    for ln, name in _fn_cache[path]:
        if ln <= line:
            best = name
        else:
            break
    return "%s:%s" % (os.path.basename(path), best)


# --------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------
def summarize(dyn):
    by_pipe, by_op, by_fn = Counter(), Counter(), Counter()
    issue = Fraction(0)
    for x, w in dyn:
        pipe, slots = pipe_of(x.op)
        issue += w
        by_pipe[pipe] += w * Fraction(slots).limit_denominator(100)
        by_op[x.op.split(".")[0]] += w
        fn = function_at(*x.chain[0]) if x.chain else x.fn[:40]
        by_fn[fn] += w
    by_pipe["issue"] = issue
    return by_pipe, by_op, by_fn


def roofline(by_pipe, machine, measured=None, clmad_rate=None):
    rows = []
    sms, ghz = machine["sms"], machine["clock_ghz"]
    per_sec = sms * ghz * 1e9
    t_meas = per_sec / (measured * 1e9) if measured else None
    for pipe, spec in machine["pipes"].items():
        rates = [(None, spec["rate"])]
        if pipe == "clmad":
            rates = list(spec.get("alternatives", {}).items()) or [(None, spec["rate"])]
            if clmad_rate:
                rates = [("--clmad-rate", clmad_rate)]
        work = float(by_pipe.get(pipe, 0))
        for label, rate in rates:
            t = work / rate
            ceiling = per_sec / t / 1e9 if t else float("inf")
            rows.append({"pipe": pipe, "rate": rate, "rateLabel": label, "work": work,
                         "smClocksPerUpdate": t, "ceilingBps": ceiling,
                         "utilisation": (t / t_meas) if t_meas else None})
    return rows, t_meas


def expected_clmad(defs):
    """CLMADs per update the arithmetic implies, for the knob combinations this
    tool knows how to count; None otherwise.  A mismatch with the dynamic count
    means a loop or branch was priced wrong."""
    d = {k: int(v) if re.fullmatch(r"-?\d+", str(v)) else v for k, v in defs.items()}
    if not d.get("ECC_WALK_TABLE") or not d.get("ECC_PACKED_CLMAD"):
        return None
    if d.get("ECC_PACKED_CLMAD_SQUARE"):
        return None
    B = d["ECC_BATCH"]
    chains = d.get("ECC_PACKED_CHAINS", 1)
    per_product = 10 if (d.get("ECC_PACKED_TOP_CLMAD") or d.get("ECC_PACKED_KARAT3")) else 6
    # 5 products per slot less 3 at each chain's first slot; 8 per inversion
    # (Itoh-Tsujii 1,2,4,...,130), or 16 when each is the two-product ONB multiply.
    inv = 16 if d.get("ECC_PACKED_ONB_INV") else 8
    products = Fraction(5 * B - 3 * chains, B) + Fraction(inv * chains, B)
    clmad = products * per_product
    if not d.get("ECC_PACKED_ALU_SQUARE"):
        clmad += 5                      # the polynomial squaring of lambda, one per update
    if not d.get("ECC_PACKED_ALU_SQR"):
        clmad += Fraction(5 * 4 * chains, B)   # five ONB squarings per inversion, four CLMADs each
    return float(clmad)


def print_report(args, defs, arch, regs, spills, dyn, m, machine):
    by_pipe, by_op, by_fn = summarize(dyn)
    B = int(defs.get("ECC_BATCH", 16))
    print("build: %s%s  arch %s  batch %d  threads %s  regs %s  spills %s" % (
        args.target, (" KNOBS=\"%s\"" % args.knobs) if args.knobs else "", arch, B,
        defs.get("ECC_THREADS"), regs, spills))
    if m.unknown:
        freq = defaultdict(float)
        for x, w in dyn:
            freq[(x.fn, x.addr)] += w
        hot = [(f, a, why, what, freq[(f, a)]) for f, a, why, what in m.unknown]
        cold = [u for u in hot if u[4] < 1e-4]
        hot = [u for u in hot if u[4] >= 1e-4]
        print("\nbranches no rule covers: %d reached less than once per 10^4 updates (inert), %d reached more often"
              % (len(cold), len(hot)))
        for fname, addr, why, what, f in hot:
            print("  %s+%#x  %.4f/update  %s  %s" % (fname[-40:], addr, f, why, what))
    exp = expected_clmad(defs)
    got = float(by_pipe.get("clmad", 0))
    if exp is not None:
        print("\nself-check: dynamic CLMAD %.3f per update, the arithmetic implies %.3f -- %s" % (
            got, exp, "OK" if abs(got - exp) < 0.01 else "MISMATCH: a loop or branch is priced wrong"))
    print("\ndynamic lane-instructions per scalar update, by pipe:")
    for pipe in ("issue", "alu", "fma", "clmad", "xu", "lsu", "const", "uniform", "control", "nop"):
        if pipe in by_pipe:
            print("  %-8s %8.2f" % (pipe, float(by_pipe[pipe])))
    if args.ops:
        print("\nby opcode:")
        for op, w in sorted(by_op.items(), key=lambda kv: -kv[1])[:args.ops]:
            print("  %-10s %8.2f" % (op, float(w)))
    if args.functions:
        print("\nby source function (innermost inlined frame):")
        for fn, w in sorted(by_fn.items(), key=lambda kv: -kv[1])[:args.functions]:
            print("  %-48s %8.2f" % (fn, float(w)))
    rows, t_meas = roofline(by_pipe, machine, args.measured, args.clmad_rate)
    print("\nroofline on %s, %d SMs at %.3f GHz (%s):" % (machine["name"], machine["sms"],
                                                        machine["clock_ghz"], machine["clock_source"]))
    hdr = "  %-7s %9s %10s %12s %12s" % ("pipe", "work/upd", "lanes/clk", "SM-clk/upd", "ceiling B/s")
    if t_meas:
        hdr += " %11s" % "busy"
    print(hdr)
    for r in rows:
        line = "  %-7s %9.1f %10.2f %12.2f %12.2f" % (r["pipe"], r["work"], r["rate"], r["smClocksPerUpdate"],
                                                     r["ceilingBps"])
        if t_meas:
            line += " %10.1f%%" % (100 * r["utilisation"])
        if r["rateLabel"]:
            line += "   " + r["rateLabel"]
        print(line)
    if t_meas:
        print("\nmeasured %.3f B/s = %.2f SM-clocks per update." % (args.measured, t_meas))
        top = max(rows, key=lambda r: r["utilisation"])
        print("busiest pipe: %s at %.1f%% of %.2f lanes/SM-clk%s." % (
            top["pipe"], 100 * top["utilisation"], top["rate"],
            (" (" + top["rateLabel"] + ")") if top["rateLabel"] else ""))
    return rows, by_pipe, by_op, by_fn


def analyse(target="gpu-rtx-pro6000-20b", knobs="", gpu="rtx-pro-6000", arch=None, steps=1024, keep=None,
            make_vars=None):
    """Build one knob set and return its dynamic per-update counts and model."""
    machine = MACHINES[gpu]
    if target == "gpu-rtx-pro6000-20b" and knobs:
        target = "gpu-preset"
    arch = arch or machine["arch"]
    arch_flag = None
    if arch != "sm_120":
        num = arch.replace("sm_", "").rstrip("a")
        arch_flag = "-gencode arch=compute_%s,code=%s" % (num, arch)
    defs, _ = make_defines(target, knobs, arch_flag, make_vars)
    work = keep or tempfile.mkdtemp()
    os.makedirs(work, exist_ok=True)
    text, regs, spills = build(defs, arch, work)
    if keep:
        with open(os.path.join(work, "walk.gi.sass"), "w") as f:
            f.write(text)
    funcs, order = parse(text, KERNEL)
    if KERNEL not in funcs:
        sys.exit("walk kernel not found in the disassembly")
    B = int(defs["ECC_BATCH"])
    threads = int(defs.get("ECC_THREADS", 256))
    dyn, m = dynamic_counts(funcs, KERNEL, B, steps, threads, 12183, int(defs.get("ECC_UNROLL_SLOTS", 1)))
    by_pipe, by_op, by_fn = summarize(dyn)
    return {"target": target, "knobs": knobs, "gpu": gpu, "arch": arch, "defs": defs, "regs": regs,
            "spills": spills, "dyn": dyn, "model": m, "by_pipe": by_pipe, "by_op": by_op, "by_fn": by_fn,
            "expectedClmad": expected_clmad(defs)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", default="gpu-rtx-pro6000-20b",
                    help="Makefile target whose knob set to analyse (gpu-rtx-pro6000-20b, gpu-preset, "
                         "gpu-b200-19b, gpu-rtx-pro6000-chains2, ...)")
    ap.add_argument("--knobs", default="", help="KNOBS= overrides, as for `make gpu-preset`")
    ap.add_argument("--make-var", action="append", default=[], help="extra make variable, e.g. CHAINS2_THREADS=384")
    ap.add_argument("--gpu", default="rtx-pro-6000", choices=sorted(MACHINES))
    ap.add_argument("--arch", help="override the compile arch (default: the machine's)")
    ap.add_argument("--measured", type=float, help="measured B complete scalar updates/s, to report busy %%")
    ap.add_argument("--clmad-rate", type=float, help="price CLMAD at this lanes/SM-clk only")
    ap.add_argument("--steps", type=int, default=1024, help="steps per launch (--bench --steps)")
    ap.add_argument("--ops", type=int, default=0, help="print the top N opcodes")
    ap.add_argument("--functions", type=int, default=0, help="print the top N source functions")
    ap.add_argument("--branches", action="store_true", help="print every loop and branch with its probability")
    ap.add_argument("--json", help="write the report here")
    ap.add_argument("--keep", help="keep the disassembly in this directory")
    a = ap.parse_args()
    machine = MACHINES[a.gpu]
    r = analyse(a.target, a.knobs, a.gpu, a.arch, a.steps, a.keep, a.make_var)
    a.target = r["target"]
    defs, m, dyn = r["defs"], r["model"], r["dyn"]
    if a.branches:
        print("loops (function tail, target, back-edge, construct, trips per entry):")
        for f, s0, e0, pat, n in m.loops:
            print("  %-24s %#7x %#7x  %s  %s" % (f[-24:], s0, e0, pat, n))
        print("forward branches (function tail, address, P(taken), regions taken side enters, other side):")
        for f, addr, p, rt, rf in m.branches:
            print("  %-24s %#7x  %.4f  %s | %s" % (f[-24:], addr, p, rt, rf))
        print()
    rows, by_pipe, by_op, by_fn = print_report(a, defs, r["arch"], r["regs"], r["spills"], dyn, m, machine)
    if a.json:
        out = {"target": a.target, "knobs": a.knobs, "gpu": a.gpu, "arch": r["arch"], "defines": defs,
               "registers": r["regs"], "spills": r["spills"], "measuredBps": a.measured,
               "perUpdate": {k: float(v) for k, v in by_pipe.items()},
               "expectedClmad": r["expectedClmad"],
               "byOpcode": {k: float(v) for k, v in by_op.items()},
               "byFunction": {k: float(v) for k, v in by_fn.items()},
               "roofline": rows, "unknownBranches": m.unknown,
               "loops": [(f[-30:], hex(s), hex(e), p, n) for f, s, e, p, n in m.loops]}
        with open(a.json, "w") as f:
            json.dump(out, f, indent=1, default=str)


if __name__ == "__main__":
    main()
