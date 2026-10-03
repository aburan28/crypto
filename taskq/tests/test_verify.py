"""The built-in verifier against published curve constants and brute force."""
import json
import subprocess
import sys

import pytest

from taskq.verify import BinaryCurve, PrimeCurve, scalar_mul, verify_certificate

# SEC 2 v2: secp256k1 and sect163k1 (NIST K-163). n*G == O exercises every
# branch of the group law hundreds of times; a wrong formula cannot pass it.
SECP256K1 = dict(
    p=2**256 - 2**32 - 977, a=0, b=7,
    G=(0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
       0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8),
    n=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141)
SECT163K1 = dict(
    m=163, modulus=[163, 7, 6, 3, 0], a=1, b=1,
    G=(0x02FE13C0537BBC11ACAA07D793DE4E6D5E5C94EEE8,
       0x0289070FB05D38FF58321F2E800536D538CCDAA3D9),
    n=0x04000000000000000000020108A2E0CC0D99F8A5EF)


def test_secp256k1_order():
    c = SECP256K1
    E = PrimeCurve(c["p"], c["a"], c["b"])
    assert E.on_curve(c["G"])
    assert scalar_mul(E, c["n"], c["G"]) is None
    assert scalar_mul(E, c["n"] - 1, c["G"]) == (c["G"][0], c["p"] - c["G"][1])


def test_sect163k1_order():
    c = SECT163K1
    E = BinaryCurve(c["m"], c["modulus"], c["a"], c["b"])
    assert E.on_curve(c["G"])
    assert scalar_mul(E, c["n"], c["G"]) is None
    x, y = c["G"]
    assert scalar_mul(E, c["n"] - 1, c["G"]) == (x, x ^ y)


def _brute_points(E, size):
    return [(x, y) for x in range(size) for y in range(size) if E.on_curve((x, y))]


@pytest.mark.parametrize("E, size", [
    (PrimeCurve(97, 2, 3), 97),
    (BinaryCurve(5, [5, 2, 0], 1, 1), 32),
    (BinaryCurve(5, [5, 2, 0], 0, 3), 32),
])
def test_toy_group_law_matches_enumeration(E, size):
    pts = _brute_points(E, size)
    N = len(pts) + 1
    for P in pts:
        assert scalar_mul(E, N, P) is None          # Lagrange
        acc = None
        for k in range(1, 12):                        # double-and-add == repeated add
            acc = E.add(acc, P)
            assert acc == scalar_mul(E, k, P)
            assert E.on_curve(acc)


def _cert(k, **over):
    c = SECP256K1
    Q = scalar_mul(PrimeCurve(c["p"], 0, 7), 123456789, c["G"])
    cert = {"kind": "discrete_log",
            "curve": {"field": "prime", "p": hex(c["p"]), "a": 0, "b": "7"},
            "statement": {"P": [hex(v) for v in c["G"]], "Q": [str(v) for v in Q],
                          "k": k, "n": hex(c["n"])}}
    cert.update(over)
    return cert


def test_certificates():
    assert verify_certificate(_cert(123456789))["status"] == "verified"
    assert verify_certificate(_cert(123456789 + SECP256K1["n"]))["status"] == "verified"
    assert verify_certificate(_cert(123456788))["status"] == "refuted"
    off = _cert(1)
    off["statement"]["Q"] = ["1", "1"]
    assert "not on the curve" in verify_certificate(off)["detail"]
    assert verify_certificate({"kind": "none"})["status"] == "no_claim"
    assert verify_certificate({"kind": "decomposition"})["status"] == "error"
    assert verify_certificate(_cert("zz"))["status"] == "error"


def test_binary_certificate_and_cli(tmp_path):
    c = SECT163K1
    E = BinaryCurve(c["m"], c["modulus"], c["a"], c["b"])
    Q = scalar_mul(E, 987654321, c["G"])
    cert = {"kind": "discrete_log",
            "curve": {"field": "binary", "m": 163, "modulus": c["modulus"], "a": 1, "b": 1},
            "statement": {"P": [hex(v) for v in c["G"]], "Q": [hex(v) for v in Q],
                          "k": "987654321"}}
    f = tmp_path / "certificate.json"
    f.write_text(json.dumps(cert))
    p = subprocess.run([sys.executable, "-m", "taskq.verify", str(f)],
                       capture_output=True, text=True)
    assert p.returncode == 0 and json.loads(p.stdout)["status"] == "verified"
    cert["statement"]["k"] = "987654322"
    f.write_text(json.dumps(cert))
    assert subprocess.run([sys.executable, "-m", "taskq.verify", str(f)]).returncode == 1
