"""Independent integer replay for the frozen P-521 Rust rows."""
from pathlib import Path
import gzip

P = (1 << 521) - 1
path = Path(__file__).with_name("p521_rows.txt.gz")
with gzip.open(path, "rt", encoding="ascii") as stream:
    for index, line in enumerate(stream):
        words = line.split()
        assert len(words) == 8, (index, len(words))
        a, b, x, y, add, sub, mul, inv = (
            int.from_bytes(bytes.fromhex(word), "little") for word in words
        )
        assert x == a % P, (index, "canonical a")
        assert y == b % P, (index, "canonical b")
        assert add == (x + y) % P, (index, "add")
        assert sub == (x - y) % P, (index, "sub")
        assert mul == (x * y) % P, (index, "mul")
        assert inv == pow(x, P - 2, P), (index, "invert")
assert index + 1 == 1000, index + 1
print("PASS: 1000 P-521 canonicalization/add/sub/mul/invert rows")
