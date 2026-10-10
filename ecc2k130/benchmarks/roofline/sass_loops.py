#!/usr/bin/env python3
"""What each pipes.cu stream's timed loop actually issues, from its SASS.

pipes.cu divides lane-instructions by SM-clocks, and the instructions it counts
are the ones its source asks for; this reads the loop body ptxas emitted for
every stream<OP> kernel so a stream the compiler rewrote cannot pass for the
pipe it was written to load.  Each of these happened while pipes.cu was being
written: a Karatsuba middle product hoisted out of the loop (32 CLMADs a round
instead of 96), FFMAs on loop-invariant data moved to the uniform datapath,
and a 32-bit mad.hi addend costing a MOV per multiply to form sm_120's 64-bit
IMAD.HI addend.

    nvcc -O3 -std=c++17 -arch=sm_120 -cubin pipes.cu -o pipes.cubin
    nvdisasm -c pipes.cubin > pipes.sass
    ./sass_loops.py pipes.sass            # one JSON line per stream

The loop is the body between the target of the one backward branch and that
branch; the per-round counts are that body divided by the 16 chains.
"""
import collections
import json
import re
import sys

CHAINS = 16


def loops(text):
    funcs = re.split(r"\n\s*\.text\.(\S+):", text)
    for i in range(1, len(funcs), 2):
        m = re.search(r"streamILi(\d+)E", funcs[i])
        if not m:
            continue
        seq, labels = [], {}
        for line in funcs[i + 1].split("\n"):
            lm = re.match(r"^\s*(\.L_x_\d+):", line)
            if lm:
                labels[lm.group(1)] = len(seq)
                continue
            im = re.search(r"/\*[0-9a-f]{4,}\*/\s+(@!?U?P\w+\s+)?([A-Z][A-Z0-9_.]*)\s*([^;]*);", line)
            if im:
                seq.append((im.group(2), im.group(3)))
        loop = None
        for k, (op, args) in enumerate(seq):
            t = re.search(r"(\.L_x_\d+)", args) if op.startswith("BRA") else None
            if t and t.group(1) in labels and labels[t.group(1)] < k:
                loop = (labels[t.group(1)], k)
        if loop is None:
            yield int(m.group(1)), None
            continue
        yield int(m.group(1)), collections.Counter(op for op, _ in seq[loop[0]:loop[1] + 1])


def main():
    text = open(sys.argv[1]).read()
    for op, mix in sorted(loops(text)):
        if mix is None:
            print(json.dumps({"op": op, "loop": None}))
            continue
        print(json.dumps({"op": op, "loopInstructions": sum(mix.values()),
                          "perChainRound": {k: round(v / CHAINS, 3) for k, v in mix.most_common()}}))


if __name__ == "__main__":
    main()
