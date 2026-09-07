#!/usr/bin/env python3
"""Differentially test the inline PTX in fp256.cuh without a GPU.

`FP_PTX=1` replaces seven multiprecision primitives with hand-written carry
chains.  Every one of those asm blocks is guarded on `__CUDA_ARCH__`, so the
host never executes them: `./bench selftest` built with `-DFP_PTX=1` is the
real check, and it needs a CUDA device.

This script covers the part of that check which does not need hardware.  It
parses the actual `asm(...)` blocks out of `fp256.cuh` -- the instruction
template and both operand lists, not a transcription of them -- binds `%N`
to the C expressions the constraint lists name, interprets the instructions
against the PTX extended-precision semantics, and compares the result with
the portable C++ branch of the same function, reimplemented here from the
`#else` side.

What that catches: operand misnumbering, a carry chain broken in the middle,
a `.cc` dropped or added, wrong limb offsets, an output read before it is
written, and bound violations in the surrounding algorithm.  Those are the
failure modes of hand-written asm and they are the ones this code has.

What it does NOT catch: anything about how the compiler allocates registers
to those operands, actual ptxas or SASS behaviour, or a mistake in the PTX
semantics encoded below.  A pass here is not a substitute for running
`./bench selftest` on a device.  It is what can be checked without one.

PTX semantics used (PTX ISA, "Extended-Precision Integer Arithmetic"):

    add.cc.u32   d, a, b      d = a + b,              CC.CF = carry out
    addc.cc.u32  d, a, b      d = a + b + CC.CF,      CC.CF = carry out
    addc.u32     d, a, b      d = a + b + CC.CF
    sub.cc.u32   d, a, b      d = a - b,              CC.CF = borrow out
    subc.cc.u32  d, a, b      d = a - b - CC.CF,      CC.CF = borrow out
    subc.u32     d, a, b      d = a - b - CC.CF
    mul.lo.u32   d, a, b      d = lo32(a * b)
    mad.lo.cc.u32  d, a, b, c   d = lo32(a*b) + c,           CC.CF = carry out
    madc.lo.cc.u32 d, a, b, c   d = lo32(a*b) + c + CC.CF,   CC.CF = carry out
    mad.hi.cc.u32  d, a, b, c   d = hi32(a*b) + c,           CC.CF = carry out
    madc.hi.cc.u32 d, a, b, c   d = hi32(a*b) + c + CC.CF,   CC.CF = carry out
    madc.hi.u32    d, a, b, c   d = hi32(a*b) + c + CC.CF

Usage:  ./ptx_asm_check.py [--trials N] [--seed S]
"""

import argparse
import os
import random
import re
import sys

M32 = 0xFFFFFFFF
HERE = os.path.dirname(os.path.abspath(__file__))

# --------------------------------------------------------------------------
# 1. Pull the asm blocks out of the header
# --------------------------------------------------------------------------


def _balanced(text, open_at):
    """Index just past the ')' matching the '(' at open_at."""
    depth, i, in_str = 0, open_at, False
    while i < len(text):
        ch = text[i]
        if in_str:
            if ch == "\\":
                i += 2
                continue
            if ch == '"':
                in_str = False
        elif ch == '"':
            in_str = True
        elif ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                return i + 1
        i += 1
    raise ValueError("unbalanced asm(")


def _split_sections(arg):
    """Split an asm argument on the top-level ':' separators."""
    out, depth, in_str, cur = [], 0, False, []
    i = 0
    while i < len(arg):
        ch = arg[i]
        if in_str:
            cur.append(ch)
            if ch == "\\":
                cur.append(arg[i + 1])
                i += 2
                continue
            if ch == '"':
                in_str = False
        elif ch == '"':
            in_str = True
            cur.append(ch)
        elif ch == "(":
            depth += 1
            cur.append(ch)
        elif ch == ")":
            depth -= 1
            cur.append(ch)
        elif ch == ":" and depth == 0:
            out.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
        i += 1
    out.append("".join(cur))
    return out


_STR = re.compile(r'"((?:[^"\\]|\\.)*)"')
_OPERAND = re.compile(r'"([^"]*)"\s*\(([^()]*(?:\([^()]*\)[^()]*)*)\)')


def _template(section):
    body = "".join(_STR.findall(section))
    return body.replace("\\n", "\n").replace("\\t", "\t")


def _operands(section):
    return [(c, e.strip()) for c, e in _OPERAND.findall(section)]


def parse_header(path):
    """-> {function name: [ {template, outputs, inputs}, ... ]} in file order."""
    src = open(path).read()
    # Strip the portable branches so only the FP_PTX side is parsed.  The
    # guards in this header are always  #if FP_PTX ... #else ... #endif.
    kept, depth_skip = [], False
    for line in src.split("\n"):
        s = line.strip()
        if s.startswith("#if FP_PTX"):
            kept.append(line)
            continue
        if s == "#else":
            depth_skip = True
            continue
        if s == "#endif":
            depth_skip = False
            continue
        if not depth_skip:
            kept.append(line)
    src = "\n".join(kept)

    funcs, order = {}, []
    fn = re.compile(r"^FP_HD\s+[\w ]+?\s+(\w+)\s*\(", re.M)
    marks = [(m.start(), m.group(1)) for m in fn.finditer(src)]
    for i, (pos, name) in enumerate(marks):
        end = marks[i + 1][0] if i + 1 < len(marks) else len(src)
        body = src[pos:end]
        blocks = []
        for m in re.finditer(r"\basm\s*\(", body):
            arg = body[m.end() - 1 : _balanced(body, m.end() - 1)][1:-1]
            secs = _split_sections(arg)
            if len(secs) < 3:
                raise ValueError(f"{name}: asm block with no input section")
            blocks.append(
                {
                    "template": _template(secs[0]),
                    "outputs": _operands(secs[1]),
                    "inputs": _operands(secs[2]),
                }
            )
        if blocks:
            funcs[name] = blocks
            order.append(name)
    return funcs, order


# --------------------------------------------------------------------------
# 2. Interpret them
# --------------------------------------------------------------------------

_LVALUE = re.compile(r"^(\w+)(?:\[(\d+)\])?$")


class AsmError(Exception):
    pass


class Slot:
    """One inline-asm operand: its own register, bound to a C lvalue."""

    __slots__ = ("expr", "constraint", "value", "defined", "is_out")

    def __init__(self, constraint, expr, env):
        self.expr, self.constraint = expr, constraint
        self.is_out = constraint[0] in "+=&" or constraint.startswith("=")
        writeonly = constraint.startswith("=")
        if writeonly:
            self.value, self.defined = None, False
        else:
            self.value, self.defined = read_lvalue(env, expr), True


def read_lvalue(env, expr):
    m = _LVALUE.match(expr)
    if not m:
        raise AsmError(f"cannot parse operand expression {expr!r}")
    name, idx = m.group(1), m.group(2)
    if name not in env:
        raise AsmError(f"unknown variable {name!r}")
    return env[name][int(idx)] if idx is not None else env[name]


def write_lvalue(env, expr, val):
    m = _LVALUE.match(expr)
    name, idx = m.group(1), m.group(2)
    if idx is not None:
        env[name][int(idx)] = val
    else:
        env[name] = val


def run_block(block, env, where):
    slots = [Slot(c, e, env) for c, e in block["outputs"]]
    n_out = len(slots)
    slots += [Slot(c, e, env) for c, e in block["inputs"]]

    def src(tok):
        tok = tok.strip()
        if tok.startswith("%"):
            s = slots[int(tok[1:])]
            if not s.defined:
                raise AsmError(
                    f"{where}: {tok} ({s.expr}, constraint {s.constraint!r}) "
                    f"is read before it is written"
                )
            return s.value
        return int(tok, 0) & M32

    def dst(tok, val):
        tok = tok.strip()
        if not tok.startswith("%"):
            raise AsmError(f"{where}: destination {tok!r} is not an operand")
        s = slots[int(tok[1:])]
        s.value, s.defined = val & M32, True

    cf = None  # condition-code carry/borrow; None until a .cc sets it
    for raw in block["template"].split(";"):
        ins = raw.strip()
        if not ins:
            continue
        op, _, rest = ins.partition(" ")
        args = [a for a in rest.split(",")]
        if len(args) < 3:
            raise AsmError(f"{where}: cannot parse {ins!r}")

        def carry_in():
            if cf is None:
                raise AsmError(f"{where}: {op} uses CC.CF before any .cc set it")
            return cf

        if op in ("add.cc.u32", "addc.cc.u32", "addc.u32"):
            t = src(args[1]) + src(args[2])
            if op != "add.cc.u32":
                t += carry_in()
            if op.endswith(".cc.u32") or op == "add.cc.u32":
                cf = t >> 32
            dst(args[0], t)
        elif op in ("sub.cc.u32", "subc.cc.u32", "subc.u32"):
            t = src(args[1]) - src(args[2])
            if op != "sub.cc.u32":
                t -= carry_in()
            if op.endswith(".cc.u32") or op == "sub.cc.u32":
                cf = 1 if t < 0 else 0
            dst(args[0], t)
        elif op == "mul.lo.u32":
            dst(args[0], (src(args[1]) * src(args[2])) & M32)
        elif op in (
            "mad.lo.cc.u32",
            "madc.lo.cc.u32",
            "mad.hi.cc.u32",
            "madc.hi.cc.u32",
            "madc.hi.u32",
        ):
            if len(args) != 4:
                raise AsmError(f"{where}: {op} needs four operands: {ins!r}")
            prod = src(args[1]) * src(args[2])
            part = (prod >> 32) if ".hi" in op else (prod & M32)
            t = part + src(args[3])
            if op.startswith("madc"):
                t += carry_in()
            if op.endswith(".cc.u32"):
                cf = t >> 32
            dst(args[0], t)
        else:
            raise AsmError(f"{where}: unmodelled instruction {op!r}")

    for s in slots[:n_out]:
        if not s.defined:
            raise AsmError(f"{where}: output {s.expr} was never written")
        write_lvalue(env, s.expr, s.value)


# --------------------------------------------------------------------------
# 3. The portable branch, transcribed from the #else sides
# --------------------------------------------------------------------------


def ref_mp_add(r, a, b):
    c = 0
    for j in range(8):
        c += a[j] + b[j]
        r[j] = c & M32
        c >>= 32
    return c


def ref_mp_sub(r, a, b):
    bw = 0
    for j in range(8):
        d = a[j] - b[j] - bw
        r[j] = d & M32
        bw = 1 if d < 0 else 0
    return bw


def ref_mp_mac_row(t, a, b):
    c = 0
    for j in range(8):
        c += t[j] + a[j] * b
        t[j] = c & M32
        c >>= 32
    c += t[8]
    t[8] = c & M32
    t[9] = (t[9] + (c >> 32)) & M32


def ref_mp_mul_row0(t, a, b):
    c = 0
    for j in range(8):
        c += a[j] * b
        t[j] = c & M32
        c >>= 32
    t[8] = c & M32


def ref_mp_mac_row9(t, a, b):
    c = 0
    for j in range(8):
        c += t[j] + a[j] * b
        t[j] = c & M32
        c >>= 32
    t[8] = c & M32


def ref_mp_add_shift32(t, a):
    c = 0
    for j in range(8):
        c += t[j + 1] + a[j]
        t[j + 1] = c & M32
        c >>= 32
    t[9] = (t[9] + c) & M32


def ref_mp_add_small(t, a0, a1, a2):
    c = t[0] + a0
    t[0] = c & M32
    c >>= 32
    c += t[1] + a1
    t[1] = c & M32
    c >>= 32
    c += t[2] + a2
    t[2] = c & M32
    c >>= 32
    for j in range(3, 8):
        c += t[j]
        t[j] = c & M32
        c >>= 32
    return c


# --------------------------------------------------------------------------
# 4. Cases: how to build an environment and what the reference should give
# --------------------------------------------------------------------------

P = (1 << 256) - (1 << 32) - 977


def limbs(n, k=8):
    return [(n >> (32 * i)) & M32 for i in range(k)]


def rnd(rng, k=8):
    return [rng.getrandbits(32) for _ in range(k)]


def edge_words(rng):
    """A word chosen to sit on a boundary more often than chance allows."""
    return rng.choice(
        [0, 1, M32, M32 - 1, 0x80000000, 0x7FFFFFFF, rng.getrandbits(32)]
    )


def edge_vec(rng, k=8):
    pick = rng.random()
    if pick < 0.15:
        return [0] * k
    if pick < 0.30:
        return [M32] * k
    if pick < 0.45:
        return limbs(P, k) if k == 8 else [M32] * k
    if pick < 0.60:
        return limbs(P - 1, k) if k == 8 else [0] * k
    if pick < 0.80:
        return [edge_words(rng) for _ in range(k)]
    return rnd(rng, k)


def case_mp_add(rng):
    a, b = edge_vec(rng), edge_vec(rng)
    env = {"r": [None] * 8, "a": list(a), "b": list(b), "c": None}
    ref = {"r": [0] * 8}
    ref["c"] = ref_mp_add(ref["r"], a, b)
    return env, ref, ("r", "c")


def case_mp_sub(rng):
    a, b = edge_vec(rng), edge_vec(rng)
    env = {"r": [None] * 8, "a": list(a), "b": list(b), "bw": None}
    ref = {"r": [0] * 8}
    ref["bw"] = ref_mp_sub(ref["r"], a, b)
    # the asm returns `bw & 1`; subc.u32 d,0,0 gives 0xFFFFFFFF on borrow
    return env, ref, ("r", "bw&1")


def case_mp_mac_row(rng):
    # mont_mul's invariant: t < 2^289 at row entry, i.e. t[9] == 0 and the
    # row cannot overflow the 10-limb window.
    t = edge_vec(rng, 8) + [rng.choice([0, 1, rng.getrandbits(16)]), 0]
    a, b = edge_vec(rng), edge_words(rng)
    env = {"t": list(t), "a": list(a), "b": b}
    ref = {"t": list(t)}
    ref_mp_mac_row(ref["t"], a, b)
    return env, ref, ("t",)


def case_mp_mul_row0(rng):
    a, b = edge_vec(rng), edge_words(rng)
    env = {"t": [None] * 9, "a": list(a), "b": b}
    ref = {"t": [0] * 9}
    ref_mp_mul_row0(ref["t"], a, b)
    return env, ref, ("t",)


def case_mp_mac_row9(rng):
    t = edge_vec(rng, 8) + [None]
    a, b = edge_vec(rng), edge_words(rng)
    env = {"t": list(t), "a": list(a), "b": b}
    ref = {"t": list(t[:8]) + [0]}
    ref_mp_mac_row9(ref["t"], a, b)
    return env, ref, ("t",)


def case_mp_add_shift32(rng):
    t = edge_vec(rng, 8) + [rng.choice([0, 1, rng.getrandbits(11)]), 0]
    a = edge_vec(rng)
    env = {"t": list(t), "a": list(a)}
    ref = {"t": list(t)}
    ref_mp_add_shift32(ref["t"], a)
    return env, ref, ("t",)


def case_mp_add_small(rng):
    t = edge_vec(rng)
    a0, a1, a2 = edge_words(rng), edge_words(rng), rng.choice([0, 1, 2])
    env = {"t": list(t), "a0": a0, "a1": a1, "a2": a2, "c": None}
    ref = {"t": list(t)}
    ref["c"] = ref_mp_add_small(ref["t"], a0, a1, a2)
    return env, ref, ("t", "c")


CASES = {
    "mp_add": case_mp_add,
    "mp_sub": case_mp_sub,
    "mp_mac_row": case_mp_mac_row,
    "mp_mul_row0": case_mp_mul_row0,
    "mp_mac_row9": case_mp_mac_row9,
    "mp_add_shift32": case_mp_add_shift32,
    "mp_add_small": case_mp_add_small,
}


# --------------------------------------------------------------------------
# 5. Drive
# --------------------------------------------------------------------------


def compare(env, ref, checks, name, trial):
    for what in checks:
        if what.endswith("&1"):
            key = what[:-2]
            got, want = env[key] & 1, ref[key] & 1
        else:
            key = what
            got, want = env[key], ref[key]
        if isinstance(want, list):
            got = [g for g in got[: len(want)]]
            want = list(want)
        if got != want:
            return (
                f"{name}: trial {trial}: {what} disagrees\n"
                f"      asm  {got}\n"
                f"      ref  {want}"
            )
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--trials", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=20260907)
    ap.add_argument("--header", default=os.path.join(HERE, "fp256.cuh"))
    args = ap.parse_args()

    funcs, order = parse_header(args.header)
    print(f"parsed {sum(len(v) for v in funcs.values())} asm blocks "
          f"in {len(funcs)} functions from {os.path.basename(args.header)}")

    missing = set(funcs) - set(CASES)
    if missing:
        print(f"  FAIL: no test case for {sorted(missing)} -- "
              f"an asm block was added without extending this script")
        return 1
    unused = set(CASES) - set(funcs)
    if unused:
        print(f"  FAIL: {sorted(unused)} has a test case but no asm block found")
        return 1

    rng = random.Random(args.seed)
    failures = 0
    for name in order:
        blocks = funcs[name]
        bad = None
        for trial in range(args.trials):
            env, ref, checks = CASES[name](rng)
            try:
                for i, blk in enumerate(blocks):
                    run_block(blk, env, f"{name}[block {i}]")
            except AsmError as e:
                bad = f"{name}: trial {trial}: {e}"
                break
            bad = compare(env, ref, checks, name, trial)
            if bad:
                break
        n = len(blocks)
        if bad:
            failures += 1
            print(f"  {name:<16} {n} block(s)   FAIL")
            print(f"      {bad}")
        else:
            print(f"  {name:<16} {n} block(s)   ok  ({args.trials} trials)")

    if failures:
        print(f"\n{failures} function(s) FAILED")
        return 1
    print(f"\nall {len(order)} asm functions agree with the portable path")
    print("NOTE: this checks the assembly's arithmetic, not register allocation")
    print("      or real device behaviour.  ./bench selftest -DFP_PTX=1 on a")
    print("      CUDA device is still the check that matters.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
