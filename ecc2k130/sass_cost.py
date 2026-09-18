#!/usr/bin/env python3
"""What each field routine costs in real SASS, and which pipe it spends.

path_cost.py answers this in clang PTX with every instruction counted as one.
Two corrections make the difference between a useful budget and a misleading
one:

  * ptxas fuses logic into LOP3 at a ratio close to two, so PTX shares are not
    what the machine issues.  This compiles with the shipping nvcc and counts
    the SASS.

  * clmad is a real sm_120 instruction but a narrow one.
    benchmarks/clmad-price/probe.cu measures 0.333 T clmad/s against 12.563 T
    LOP3/s on an RTX PRO 4500: one clmad occupies its unit for as long as
    37.7 LOP3 occupy the integer pipe.  Counting it as one instruction, as
    THROUGHPUT-30B.md's budget does, understates the multiplier ten-fold.

The two pipes are reported apart because they overlap: a kernel spends
max(alu, clmad) on them, not the sum.  Which one binds is not this script's
to decide, and the pure-stream rates above must not be used to decide it --
an independent LOP3 stream reaches issue rates a walk full of dependencies
cannot, so dividing by 12.563 T flatters the kernel.  Nsight Compute is the
authority: at the audited preset it reports the ALU pipe at 87.3% and the
FP64 pipe that carries clmad at 51.4%, so the walk is ALU-bound with
carryless headroom.

Read the pipe lines below as marginal exchange rates, not as headroom.  What
they support is the trade: removing a percent of the ALU column is worth
about a percent of throughput while the ALU pipe binds, and a clmad bought
with 38 logic ops pays until the carryless unit catches up.

    ./sass_cost.py                       # per-routine SASS, shipping preset
    ./sass_cost.py --update              # weighted into one scalar update
    ./sass_cost.py --rate 5.022 --sms 82 # what each pipe then costs per update

It needs nvcc 13.3 or newer, for clmad.  No pip wheel carries one; NVIDIA's apt
repository serves the .deb files and a .deb is an ar archive, so about 40 MB of
cuda-nvcc, libnvvm, cuda-crt, cuda-cudart-dev, cuda-cuobjdump and cuda-nvdisasm
unpack without root -- TOP-CLMAD.md has the six lines.  Then --nvcc and
--cuobjdump point at them.

Two things this counted wrong until it was run that way.  It matched routine
names as substrings, so "k_mul" also matched k_mulPair and k_mulPoly and the
inverse-chain row carried a 164-instruction routine that costs 347.  And it
costed the walk's Frobenius network as the single-coordinate form twice, 506,
where the pair the walk actually calls shares its mask stream and is 469.
Both are what happens to a tool nobody can run: it is trusted instead.
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))

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

# Measured by benchmarks/clmad-price/probe.cu on an RTX PRO 4500.
LOP3_RATE = 12.563e12          # lane-ops/s, the integer pipe
CLMAD_RATE = 0.333e12          # lane-ops/s, the carryless unit
CLMAD_PRICE = LOP3_RATE / CLMAD_RATE
IMAD_PRICE = 2.01

# One kernel per routine: ptxas inlines device functions regardless of
# __noinline__, so a separate entry point is the only way to see a routine's
# own SASS.  k_nop calibrates the parameter load and store every kernel here
# carries, and is subtracted from the rest.
SPIKE = r"""
#include "../include/curveparams.h"
#include "../include/packedkernels.cuh"
using namespace eccPacked131;
#define KERNEL(name, expr) \
  __global__ void name(P131 *o, const P131 *i, uint32_t *h, int *w, int k) { \
      (void)h; (void)w; (void)k; *o = (expr); }
KERNEL(k_nop, i[0])
KERNEL(k_toPoly, toPolynomial131(i[0]))
KERNEL(k_fromPoly, fromPolynomial131(i[0]))
KERNEL(k_sqr, sqr131(i[0]))
KERNEL(k_sqrPoly, squarePolynomial131(i[0]))
KERNEL(k_add, add131(i[0], i[1]))
KERNEL(k_mulPoly, mulPolynomial131(i[0], i[1]))
KERNEL(k_mul, mul131(i[0], i[1]))
KERNEL(k_inv, inv131(i[0]))
KERNEL(k_sigmaWalk, sigma131(i[0], 3 + (k & 7)))
KERNEL(k_reduce, reducePolynomial131(h))
KERNEL(k_fromPolyProduct, fromPolynomialProduct131(h))
__global__ void k_product(P131 *o, const P131 *i, uint32_t *h, int *w, int k) {
    (void)w; (void)k;
    uint32_t c[9];
    product131(i[0], i[1], c);
#pragma unroll
    for (int j = 0; j < 9; ++j) h[j] = c[j];
    (void)o;
}
__global__ void k_weight(P131 *o, const P131 *i, uint32_t *h, int *w, int k) {
    (void)o; (void)h; (void)k;
    *w = weight(i[0]);
}
__global__ void k_mulOnb(P131 *o, const P131 *i, uint32_t *h, int *w, int k) {
    (void)h; (void)w; (void)k;
    *o = mulOnb131(i[0], i[1]);
}
__global__ void k_mulPair(P131 *o, const P131 *i, uint32_t *h, int *w, int k) {
    (void)h; (void)w; (void)k;
    PolynomialPair p = mulPolynomialPair131(i[0], i[1], i[2]);
    o[0] = p.first; o[1] = p.second;
}
#if ECC_PACKED_WEIGHTED_PREFIX == 2
// The walk applies the network to both coordinates through one mask stream;
// costing the single-coordinate form twice overstates it (506 against 469).
__global__ void k_sigmaPair(P131 *o, const P131 *i, uint32_t *h, int *w, int k) {
    (void)h; (void)w;
    SigmaWalkPair131 p = sigmaWalkNetworkPair131(i[0], i[1], k & 7);
    o[0] = p.first; o[1] = p.second;
}
#endif
// The 4x4-word product alone, without the three-bit top-word correction.
// Wrong as arithmetic; it isolates what that correction costs, which at the
// shipping preset is 65 of product131's 77 instructions.  See TOP-CLMAD.md.
__global__ void k_productNoTop(P131 *o, const P131 *i, uint32_t *h, int *w, int k) {
    (void)w; (void)k; (void)o;
    uint32_t c[9];
    clmul128(c, i[0].v, i[1].v); c[8] = 0;
#pragma unroll
    for (int j = 0; j < 9; ++j) h[j] = c[j];
}
"""

ALU = ("LOP3", "LOP", "IADD3", "IADD", "SHF", "SHL", "SHR", "PRMT", "POPC", "FLO",
       "BREV", "SEL", "ISETP", "IMNMX", "MOV", "PLOP3", "P2R", "R2P", "VOTE",
       "IABS", "LEA", "XOR", "BFE", "BFI", "BMSK", "FSEL", "F2I", "I2F", "IADD32I")
MEM = ("LDG", "STG", "LDS", "STS", "LDL", "STL", "LD", "ST", "LDC", "RED", "ATOM",
       "ATOMG", "MEMBAR", "CCTL", "LDSM", "ULDC")
CTRL = ("BRA", "BSYNC", "BSSY", "EXIT", "CALL", "RET", "NOP", "BAR", "SSY", "SYNC",
        "JMP", "PBK", "BRK", "WARPSYNC", "YIELD", "S2R", "CS2R", "DEPBAR", "BMOV")


def classify(op):
    head = op.split(".")[0]
    if head.startswith("CLMAD") or head.startswith("CLMUL"):
        return "clmad"
    # Blackwell issues register moves and some adds as IMAD forms; those are
    # ordinary pipe ops, not multiplies.  Only a real multiply costs two.
    if op.startswith(("IMAD.MOV", "IMAD.IADD", "IMAD.SHL")):
        return "alu"
    if head.startswith("IMAD") or head.startswith("IMUL"):
        return "imad"
    if head in MEM:
        return "mem"
    if head in CTRL:
        return "ctrl"
    return "alu"


def build(defs, arch, nvcc, work):
    src = os.path.join(work, "sass_cost_spike.cu")
    open(src, "w").write(SPIKE)
    # The spike sits in a temp dir; point the relative includes back at the tree.
    cubin = os.path.join(work, "spike.cubin")
    cmd = [nvcc, "-O3", "-std=c++17", "-arch=" + arch, "-cubin", "-I", HERE,
           "-I", os.path.join(HERE, "include"), "-o", cubin, src]
    cmd[1:1] = ["-D%s=%s" % (k, v) for k, v in defs.items()]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=os.path.join(HERE, "src"))
    if r.returncode:
        sys.exit("nvcc failed:\n" + r.stderr[:3000])
    return cubin


def functions(cubin, cuobjdump):
    r = subprocess.run([cuobjdump, "-sass", cubin], capture_output=True, text=True)
    if r.returncode:
        sys.exit("cuobjdump failed:\n" + r.stderr[:2000])
    out, name = {}, None
    for line in r.stdout.splitlines():
        m = re.match(r"\s*Function : (\S+)", line)
        if m:
            name = m.group(1)
            out[name] = Counter()
            continue
        if name is None:
            continue
        m = re.match(r"\s*/\*[0-9a-f]{4,}\*/\s+(?:@!?U?P\d\s+)?([A-Z][A-Z0-9_.]*)", line)
        if m:
            out[name][m.group(1)] += 1
    return out


def cost(counter):
    k = Counter()
    for op, n in counter.items():
        k[classify(op)] += n
    k["aluSlots"] = k["alu"] + k["imad"] * IMAD_PRICE
    k["clmadSlots"] = k["clmad"] * CLMAD_PRICE
    return k


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arch", default="sm_120")
    ap.add_argument("--nvcc", default="nvcc")
    ap.add_argument("--cuobjdump", default="cuobjdump")
    ap.add_argument("--batch", type=int, default=PRESET["ECC_BATCH"])
    ap.add_argument("--update", action="store_true", help="weight into one scalar update")
    ap.add_argument("--rate", type=float, default=0.0, help="measured B updates/s")
    ap.add_argument("--sms", type=int, default=82, help="SMs on the measured part")
    ap.add_argument("--define", action="append", default=[], help="extra -D for the spike")
    a = ap.parse_args()

    defs = dict(PRESET, ECC_BATCH=a.batch)
    for d in a.define:
        k, _, v = d.partition("=")
        defs[k] = v or 1
    work = tempfile.mkdtemp()
    fn = functions(build(defs, a.arch, a.nvcc, work), a.cuobjdump)

    order = ["k_product", "k_productNoTop", "k_reduce", "k_fromPolyProduct", "k_mulPoly", "k_mulPair",
             "k_mulOnb", "k_mul", "k_inv", "k_sigmaWalk", "k_sigmaPair", "k_toPoly", "k_fromPoly", "k_sqr",
             "k_sqrPoly", "k_add", "k_weight"]
    named = {}
    # Itanium mangling puts the name's length directly before it: _Z5k_mulP...
    # against _Z9k_mulPolyP...  A substring test had "k_mul" matching k_mulPair
    # and k_mulPoly as well, so the inverse-chain row reported whichever of the
    # three cuobjdump listed last; and "k_mulP" as a prefix does the same, since
    # Poly and Pair both start with P.  The length is what makes it exact.
    for key in fn:
        for name in ["k_nop"] + order:
            if re.search(r"_Z%d%sP" % (len(name), re.escape(name)), key):
                named[name] = fn[key]
    overhead = cost(named["k_nop"])["aluSlots"] if "k_nop" in named else 0.0
    print("SASS per routine, %s, batch %d  (less %.0f slots of kernel overhead)"
          % (a.arch, a.batch, overhead))
    print("%-22s %6s %6s %6s %6s %10s" % ("routine", "alu", "imad", "clmad", "mem", "aluSlots"))
    got = {}
    for name in order:
        if name not in named:
            continue
        k = cost(named[name])
        k["aluSlots"] = max(0.0, k["aluSlots"] - overhead)
        got[name] = k
        print("%-22s %6d %6d %6d %6d %10.0f"
              % (name, k["alu"], k["imad"], k["clmad"], k["mem"], k["aluSlots"]))

    if not a.update:
        return
    B = a.batch
    pair = "k_sigmaPair" in got
    rows = [("polynomial product, paired", "k_mulPair", (B - 1) / B),
            ("polynomial product, single", "k_mulPoly", (B + 1) / B),
            ("normal-basis product, inverse chain", "k_mul", 8.0 / B),
            ("Frobenius network, walk (pair)" if pair else "Frobenius network, walk (2 coords)",
             "k_sigmaPair" if pair else "k_sigmaWalk", 1.0 if pair else 2.0),
            ("basis conversion out (2 coords)", "k_fromPoly", 2.0),
            ("basis conversion in (2 coords)", "k_toPoly", 2.0),
            ("polynomial squaring", "k_sqrPoly", 1.0),
            ("inversion chain", "k_inv", 1.0 / B),
            ("class weight", "k_weight", 1.0),
            ("additions", "k_add", 4.0)]
    print("\nper scalar update, batch %d" % B)
    print("%-38s %6s %10s %10s" % ("", "per", "aluSlots", "clmad"))
    alu = clm = 0.0
    parts = []
    for label, name, mult in rows:
        if name not in got:
            continue
        k = got[name]
        alu += k["aluSlots"] * mult
        clm += k["clmad"] * mult
        parts.append((label, k["aluSlots"] * mult))
        print("%-38s %6.2f %10.0f %10.1f" % (label, mult, k["aluSlots"] * mult, k["clmad"] * mult))
    print("%-38s %6s %10.0f %10.1f" % ("TOTAL of the routines above", "", alu, clm))
    print()
    for label, v in sorted(parts, key=lambda p: -p[1]):
        print("   %-44s %7.0f  %5.1f%%" % (label, v, 100 * v / alu))

    if a.rate:
        print("\n--- each pipe at %.3f B updates/s on %d SMs ---" % (a.rate, a.sms))
        aluUse = a.rate * 1e9 * alu
        clmUse = a.rate * 1e9 * clm
        print("integer pipe : %7.2f T slots/s, %5.1f%% of an independent LOP3 stream"
              % (aluUse / 1e12, 100 * aluUse / LOP3_RATE))
        print("carryless    : %7.2f T clmad/s, %5.1f%% of an independent clmad stream"
              % (clmUse / 1e12, 100 * clmUse / CLMAD_RATE))
        # The literal percentages are escaped: this string is a format string,
        # and "87.3%," was being read as a conversion, which made --rate crash
        # at the last line after printing everything above it.
        print("\nThese fractions are not headroom: an independent stream issues at a rate"
              "\na dependent walk does not reach.  Nsight Compute measures the pipes"
              "\nthemselves -- ALU 87.3%%, FP64 (clmad) 51.4%% at this preset -- so the ALU"
              "\ncolumn is what binds, and one clmad is worth %.0f logic ops until the"
              "\ncarryless unit closes that gap." % CLMAD_PRICE)


if __name__ == "__main__":
    main()
