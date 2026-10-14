"""List conditional branch instructions in the frozen ARM64 assembly."""

import re
from pathlib import Path

branch = re.compile(r"\s*(b\.[a-z]+|cbz|cbnz|tbz|tbnz)\b")

for filename in ("arm64_x25519_before_barrier.s", "arm64_x25519.s"):
    assembly = (Path(__file__).parent / filename).read_text().splitlines()
    symbols = [i for i, line in enumerate(assembly) if line.startswith("__ZN") and line.endswith(":")]
    for needle in ("Fe3mul", "6x255196x25519"):
        start = next(i for i in symbols if needle in assembly[i])
        end = next((i for i in symbols if i > start), len(assembly))
        found = [(i + 1, assembly[i].strip()) for i in range(start, end) if branch.match(assembly[i])]
        print(f"{filename} {needle}: {len(found)} conditional branches")
        for line, instruction in found:
            print(f"  {line}: {instruction}")
        if needle == "Fe3mul":
            assert not found
        elif filename == "arm64_x25519.s":
            assert len(found) == 9 and all(instruction.startswith("b.ne") for _, instruction in found)
        else:
            assert len(found) == 11
            assert sum(instruction.startswith("b.lo") for _, instruction in found) == 2
