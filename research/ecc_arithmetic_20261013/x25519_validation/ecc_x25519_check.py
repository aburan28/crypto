from pathlib import Path

P = 2**255 - 19
rows = (Path(__file__).parent / "x25519_rows.txt").read_text().splitlines()
assert len(rows) == 100


def reference(scalar_bytes, u_bytes):
    k = bytearray(scalar_bytes)
    k[0] &= 248
    k[31] &= 127
    k[31] |= 64
    scalar = int.from_bytes(k, "little")
    u = int.from_bytes(u_bytes, "little") & ((1 << 255) - 1)
    u %= P
    x2, z2, x3, z3, swap = 1, 0, u, 1, 0
    for t in range(254, -1, -1):
        bit = (scalar >> t) & 1
        swap ^= bit
        if swap:
            x2, x3 = x3, x2
            z2, z3 = z3, z2
        swap = bit
        a = (x2 + z2) % P
        aa = a * a % P
        b = (x2 - z2) % P
        bb = b * b % P
        e = (aa - bb) % P
        c = (x3 + z3) % P
        d = (x3 - z3) % P
        da = d * a % P
        cb = c * b % P
        x3 = (da + cb) ** 2 % P
        z3 = u * ((da - cb) ** 2) % P
        x2 = aa * bb % P
        z2 = e * (aa + 121665 * e) % P
    if swap:
        x2, x3 = x3, x2
        z2, z3 = z3, z2
    return (x2 * pow(z2, P - 2, P) % P).to_bytes(32, "little")


for index, row in enumerate(rows):
    scalar_hex, u_hex, answer_hex = row.split()
    expected = reference(bytes.fromhex(scalar_hex), bytes.fromhex(u_hex))
    assert bytes.fromhex(answer_hex) == expected, index
print(f"verified {len(rows)} X25519 outputs against independent integer ladder")
