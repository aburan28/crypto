#!/usr/bin/env python3
"""Where the packed walk's instructions go, without a GPU.

Compiles the packed walk kernel offline with clang and counts the PTX body of
every field routine it calls, then weights each by how often a scalar update
runs it.  The point is the *composition*: which parts of the update are worth
attacking, and what a throughput objective implies for each of them.

Two things this is not.  The unit is clang's PTX, not the SASS the shipping
nvcc build executes, so absolute counts here are not the "instruction visits
per scalar update" the benchmark receipts quote -- ptxas fuses logic into
LOP3 and the ratio is close to two.  Shares are the robust output; absolute
scale comes from calibrating against a measured rate.

And ALU and memory instructions are counted apart on purpose.  The
benchmarks/vector4-sigma screen replaced 56 scalar LDS with 14 LDS.128, cut
1.9% of instruction visits and measured +0.285% throughput: memory-pipe
instructions are nearly free at the margin because the integer pipe is what
saturates.  Only the ALU column predicts speed.

    ./path_cost.py                 # shipping RTX PRO 6000 preset
    ./path_cost.py --clmad 0       # the software-product path it replaced
    ./path_cost.py --rate 14.6375  # scale the composition to a measured B/s
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))

# The audited RTX PRO 6000 preset, as RTX_PRO6000_* in the Makefile.
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

SPIKE = """
#include "../include/curveparams.h"
#include "../include/packedkernels.cuh"
namespace eccPacked131 {
/* The walk inlines these; wrap them so their bodies can be counted apart. */
__device__ __noinline__ P131 n_toPoly(P131 a) { return toPolynomial131(a); }
__device__ __noinline__ P131 n_fromPoly(P131 a) { return fromPolynomial131(a); }
__device__ __noinline__ P131 n_sqr(P131 a) { return sqr131(a); }
__device__ __noinline__ P131 n_add(P131 a, P131 b) { return add131(a, b); }
__device__ __noinline__ int n_weight(P131 a) { return weight(a); }
__global__ void probe(P131 *o, const P131 *i) {
    o[0] = n_toPoly(i[0]);  o[1] = n_fromPoly(i[1]);
    o[2] = n_sqr(i[2]);     o[3] = n_add(i[3], i[4]);
    o[4].v[0] = n_weight(i[5]);
}
}
"""

# Memory, control and call opcodes: they issue on the LSU or the branch unit,
# not the integer pipe the kernel saturates.
NON_ALU = ("ld", "st", "call", "bra", "cvta", "ret", "bar", "atom", "red",
           "(", ");", "@", "prefetch", "membar")


def parse_ptx(path):
    """Instruction lists for every function *definition* in a PTX module."""
    txt = open(path).read()
    best = {}
    for m in re.finditer(r"^\s*(?:\.visible\s+|\.weak\s+)?\.(func|entry)\b", txt, re.M):
        nm = re.search(r"(_Z\w+)\s*\(", txt[m.end():m.end() + 800])
        if not nm:
            continue
        rest = txt[m.end() + nm.end():]
        brace, semi = rest.find("{"), rest.find(";")
        if brace < 0 or (0 <= semi < brace):
            continue          # a prototype, not a definition
        i = m.end() + nm.end() + brace
        j, d = i, 0
        while True:
            if txt[j] == "{":
                d += 1
            elif txt[j] == "}":
                d -= 1
                if d == 0:
                    break
            j += 1
        ops = [l.strip() for l in txt[i:j].split("\n")
               if l.strip() and not l.strip().startswith(("//", ".", "$", "{", "}"))]
        if nm.group(1) not in best or len(ops) > len(best[nm.group(1)]):
            best[nm.group(1)] = ops
    return best


def split_alu(ops):
    alu = sum(1 for o in ops if not o.split()[0].split(".")[0].startswith(NON_ALU))
    return alu, len(ops) - alu


def compile_ptx(defs, work, clang, cuda_path):
    src = os.path.join(HERE, "src", "path_cost_spike.cu")
    with open(src, "w") as f:
        f.write(SPIKE)
    out = os.path.join(work, "spike.ptx")
    cmd = [clang, "-x", "cuda", "--cuda-device-only", "--cuda-gpu-arch=sm_90",
           "--cuda-path=" + cuda_path, "-Wno-unknown-cuda-version",
           "-Wno-macro-redefined", "-O3", "-std=c++17",
           # clang reports the pip wheel's CUDA version; clmad's guard wants 13.3.
           "-U__CUDACC_VER_MAJOR__", "-U__CUDACC_VER_MINOR__",
           "-D__CUDACC_VER_MAJOR__=13", "-D__CUDACC_VER_MINOR__=3",
           "-S", "-o", out, src]
    cmd[-3:-3] = ["-D%s=%s" % (k, v) for k, v in defs.items()]
    r = subprocess.run(cmd, capture_output=True, text=True)
    os.unlink(src)
    if r.returncode:
        sys.exit("clang failed:\n" + r.stderr[:2000])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--clmad", type=int, default=1, help="native carryless multiply")
    ap.add_argument("--batch", type=int, default=PRESET["ECC_BATCH"])
    ap.add_argument("--rate", type=float, default=14.637530,
                    help="measured B complete scalar iterations/s to scale to")
    ap.add_argument("--target", type=float, default=30.0, help="objective, B/s")
    ap.add_argument("--clang", default=os.environ.get("CLANG", "clang++"))
    ap.add_argument("--cuda-path", default=os.environ.get("CUDAPATH", "/usr/local/cuda"))
    a = ap.parse_args()

    if not shutil.which(a.clang):
        sys.exit("no %s on PATH; set CLANG" % a.clang)
    defs = dict(PRESET, ECC_PACKED_CLMAD=a.clmad, ECC_BATCH=a.batch)
    B = a.batch

    work = tempfile.mkdtemp()
    try:
        fn = parse_ptx(compile_ptx(defs, work, a.clang, a.cuda_path))
    finally:
        shutil.rmtree(work, ignore_errors=True)

    def body(frag):
        hits = [k for k in fn if frag in k]
        if not hits:
            sys.exit("no PTX body matching %r" % frag)
        return fn[max(hits, key=lambda k: len(fn[k]))]

    walk = body("L4walkE")
    calls = Counter(re.findall(r"call\.uni[^;]*?,\s*(_Z\w+)", "\n".join(walk), re.S))

    def called(frag):
        hits = [k for k in calls if frag in k]
        return calls[hits[0]] if hits else 0

    # How often a *scalar update* runs each routine.  The walk keeps two
    # unrolled-once slot loops inside one step, so a batch of B updates runs
    # each slot-loop body B times and the inversion chain once.
    #   forward : slot 0 seeds the prefix, slots 1..B-1 take a paired product
    #   backward: slots 1..B-1 take a paired product, slot 0 a single one,
    #             and every slot takes the addition's single product
    #   inverse : Itoh-Tsujii, once per batch
    rows = [
        ("polynomial product, paired (2 products)", "mulPolynomialPair131",
         2 * (B - 1) / B),
        ("polynomial product, single", "mulPolynomial131", (B + 1) / B),
        ("normal-basis product, inverse chain", "L6mul131", called("L6mul131") / B),
        ("Frobenius network, walk (2 coords)", "PairShared131", 1.0),
        ("Frobenius network, inverse chain", "sigmaInvNetwork131",
         called("sigmaInvNetwork131") / B),
        ("Frobenius network, inverse chain (walk form)", "L19sigmaWalkNetwork131",
         called("L19sigmaWalkNetwork131") / B),
    ]

    print("packed walk composition, batch %d, clmad %d" % (B, a.clmad))
    print("unit: clang PTX instructions per scalar update (not SASS visits)\n")
    print("%-46s %6s %8s %8s %8s" % ("", "per", "ALU", "mem", "ALU/upd"))
    alu_tot = mem_tot = 0.0
    parts = []
    for label, frag, mult in rows:
        if mult == 0:
            continue
        alu, mem = split_alu(body(frag))
        alu_tot += alu * mult
        mem_tot += mem * mult
        parts.append((label, alu * mult))
        print("%-46s %6.2f %8d %8d %8.0f" % (label, mult, alu, mem, alu * mult))

    # Everything the walk inlines: basis conversions, the point addition's
    # additions and squaring, the class weight, compact-state addressing, the
    # loop and the distinguished-point test.  Counted as the kernel body less
    # nothing, because every one of them is already inside it.
    alu, mem = split_alu(walk)
    print("%-46s %6s %8d %8d %8s" % ("inlined in the kernel body", "(body)", alu, mem, "-"))
    print("   of which, measured separately:")
    for label, frag in [("toPolynomial131", "n_toPoly"), ("fromPolynomial131", "n_fromPoly"),
                        ("sqr131", "n_sqr"), ("add131", "n_add"),
                        ("class weight (5 popc)", "n_weight")]:
        ia, im = split_alu(body(frag))
        print("      %-40s %14d %8d" % (label, ia, im))

    # The inlined body is one forward-loop body, one backward-loop body and the
    # inversion's own inline code; the slot loops run B times each.
    inline_per_update = alu * (2.0 / 3.0)      # two of its three regions are per-slot
    alu_tot += inline_per_update
    parts.append(("inlined per-slot work", inline_per_update))
    print("\n%-46s %6s %8s %8s %8.0f" %
          ("per-slot share of the inlined body", "~2/3", "", "", inline_per_update))

    print("\n%-46s %30.0f ALU PTX/update" % ("TOTAL", alu_tot))
    print("%-46s %30.0f mem PTX/update (near-free at the margin)" % ("", mem_tot))
    print()
    for label, v in sorted(parts, key=lambda p: -p[1]):
        print("   %-50s %6.0f  %5.1f%%" % (label, v, 100 * v / alu_tot))

    print("\n--- what the objective costs ---")
    if a.clmad != PRESET["ECC_PACKED_CLMAD"] or a.batch != PRESET["ECC_BATCH"]:
        print("NOTE: --rate defaults to the shipping preset's measured rate;"
              " pass the rate this configuration measured.")
    print("measured %.4f B/s at this composition" % a.rate)
    need = alu_tot * a.rate / a.target
    print("%.0f B/s needs %.0f ALU PTX/update: a %.2fx cut" % (a.target, need, alu_tot / need))
    prod = sum(v for l, v in parts if "product" in l)
    print("the %.0f products/update alone are %.0f ALU PTX (%.0f%% of that budget)"
          % (2 * (B - 1) / B * 2 + (B + 1) / B + called("L6mul131") / B,
             prod, 100 * prod / need))
    print("so the ceiling with every non-product instruction removed is %.1f B/s"
          % (a.rate * alu_tot / prod))


if __name__ == "__main__":
    main()
