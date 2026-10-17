from pathlib import Path

P = 2**255 - 19
rows = (Path(__file__).parent / "field_rows.txt").read_text().splitlines()
assert len(rows) == 1000
for index, row in enumerate(rows):
    a_hex, b_hex, add_hex, sub_hex, mul_hex, inv_hex = row.split()
    a = int.from_bytes(bytes.fromhex(a_hex), "little") & ((1 << 255) - 1)
    b = int.from_bytes(bytes.fromhex(b_hex), "little") & ((1 << 255) - 1)
    a %= P
    b %= P
    results = [
        int.from_bytes(bytes.fromhex(value), "little")
        for value in (add_hex, sub_hex, mul_hex, inv_hex)
    ]
    expected = [(a + b) % P, (a - b) % P, (a * b) % P, pow(a, P - 2, P)]
    assert results == expected, (index, results, expected)
print(f"verified {len(rows)} fixed-radix field rows against Python integers")
