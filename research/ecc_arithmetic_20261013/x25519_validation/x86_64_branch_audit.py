"""List conditional jumps in frozen x86-64 Rust arithmetic functions."""
from collections import Counter
from pathlib import Path

path = Path(__file__).with_name("x86_64_x25519.s")
lines = path.read_text().splitlines()


def function(name: str) -> list[str]:
    starts = [i for i, line in enumerate(lines) if line.endswith(":") and name in line]
    assert len(starts) == 1, (name, starts)
    start = starts[0]
    ends = [i for i in range(start + 1, len(lines)) if lines[i].startswith("\t.size\t")]
    assert ends, name
    return lines[start : ends[0]]


for name, expected in (("Fe3mul", {}), ("x255196x25519", {"jb": 1, "jne": 8})):
    body = function(name)
    jumps = Counter(
        line.strip().split()[0]
        for line in body
        if line.lstrip().startswith("j") and not line.lstrip().startswith(("jmp", "jmpq"))
    )
    assert dict(jumps) == expected, (name, jumps)
    print(f"{name}: {dict(jumps)}")
