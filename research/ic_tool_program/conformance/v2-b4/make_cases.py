#!/usr/bin/env python3
"""Conformance cases for B4: binary fields of three words and more
(`research/ic_tool_program/rounds/B4-multi-word/PROTOCOL.md`).

    python3 make_cases.py           # writes params/, cases.json and SHA256SUMS; refuses to overwrite
    python3 make_cases.py --check   # re-derives every file and compares the bytes

The documents are built with B1's generator code (`../v2/make_cases.py`),
whose arithmetic shares nothing with the Rust tool:
- each curve's order comes from the trace recurrence, and `r` and `h` are
  stated here and checked against it, with `r`'s primality exact;
- each instance's signed Frobenius group acts freely on its subgroup: the
  eigenvalue of Frobenius there has order exactly `n` (design §4);
- each generator is `[h]P` for a point hashed from a public label;
- each target is a known multiple of the generator, the multiplier hashed
  from a public label, so every case checks an exact answer.

C097's expected target is derived here by design §2's rule, with a BLAKE3
written in Python below, so that the tool's derivation is checked against
code that shares nothing with it.

`../run.py` runs every step's cases, these among them, with the `until`
rule applied.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
V2 = HERE.parent / "v2"
B3B = HERE.parent / "v2-b3b"
_spec = importlib.util.spec_from_file_location("conformance_v2_cases", V2 / "make_cases.py")
v2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(v2)

STEP = "B4"
RHO_SEED = 0x230000 + 1  # the suite's rho seed for target T01, as B1, B3 and B3b use
GATE_MODULUS = (1 << 83) | (1 << 45) | (1 << 2) | (1 << 1) | 1  # AGENTS.md §8a
LABEL = "ic-conformance-v2-b4"
GATE = "{cases}/gate-m83-T001.json"
SMOKE = "{cases}/C009-translation-of-smoke-row.json"
CHALLENGE = "{cases}/ecc2k130-challenge.json"
C080_DOC = "C080-kic-two-word-n79.json"  # B3b's, copied byte for byte (C095)
C078_DOC = "C078-kic-two-word-n67.json"  # B3b's, copied byte for byte (C100, C101)
TIMEOUT = 1800

# (case id, a, n, r, purpose). The modulus is the repository's rule for
# each degree (docs/curves/ICV1.md), and h = #K_a / r.
INSTANCES = [
    ("C088", 0, 127, 1118452171, "the smallest three-word degree"),
    ("C089", 0, 137, 191818802977, "a 37-bit subgroup at n = 137"),
    ("C090", 0, 151, 28016550571, "n = 151, a trinomial field"),
    ("C091", 1, 157, 154781804431543, "the largest subgroup, 47 bits"),
    ("C092", 1, 173, 208110697, "the smallest subgroup, under a 146-bit cofactor"),
    ("C093", 0, 179, 820651535909, "the widest three-word degree here, n = 179"),
]

# C098: past nine words. #K_0(GF(2^577)) = 4 * 2347237 * (a composite
# cofactor), so r = 2347237 names a subgroup.  Schema v2's binary degrees
# stop at 571 (its §3), inside nine words' 574, so the document is
# refused as `degree-range` before any gate is read; the width gate past
# nine words is the router's own test.
WIDE_N, WIDE_A, WIDE_R = 577, 0, 2347237

# C102: the 95-bit prime factor of #K_0(GF(2^127)), and the bases of its
# primality screen.
R95 = 38030500514909642861048567999
PRIMES_40 = [p for p in range(2, 180) if all(p % d for d in range(2, int(p ** 0.5) + 1))][:40]

# C099: sect163k1 (SEC 2), K_1 over x^163 + x^7 + x^6 + x^3 + 1, with its
# published generator, r and cofactor 2.
SECT163K1 = {
    "modulus": (1 << 163) | (1 << 7) | (1 << 6) | (1 << 3) | 1,
    "r": 0x4000000000000000000020108A2E0CC0D99F8A5EF,
    "h": 2,
    "g": (0x02FE13C0537BBC11ACAA07D793DE4E6D5E5C94EEE8, 0x0289070FB05D38FF58321F2E800536D538CCDAA3D9),
}


# -- BLAKE3, for C097 ------------------------------------------------------------------------
# The reference algorithm, with the extended output design §2 reads. It is
# checked against BLAKE3's published test vectors in selftest().

_IV = (0x6A09E667, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A, 0x510E527F, 0x9B05688C, 0x1F83D9AB, 0x5BE0CD19)
_PERMUTATION = (2, 6, 3, 10, 7, 0, 4, 13, 1, 11, 12, 5, 9, 14, 15, 8)
_CHUNK_START, _CHUNK_END, _PARENT, _ROOT = 1, 2, 4, 8


def _rotr(v: int, s: int) -> int:
    return ((v >> s) | (v << (32 - s))) & 0xFFFFFFFF


def _g(st: list[int], a: int, b: int, c: int, d: int, x: int, y: int) -> None:
    st[a] = (st[a] + st[b] + x) & 0xFFFFFFFF
    st[d] = _rotr(st[d] ^ st[a], 16)
    st[c] = (st[c] + st[d]) & 0xFFFFFFFF
    st[b] = _rotr(st[b] ^ st[c], 12)
    st[a] = (st[a] + st[b] + y) & 0xFFFFFFFF
    st[d] = _rotr(st[d] ^ st[a], 8)
    st[c] = (st[c] + st[d]) & 0xFFFFFFFF
    st[b] = _rotr(st[b] ^ st[c], 7)


def _compress(cv, words: list[int], counter: int, block_len: int, flags: int) -> list[int]:
    st = list(cv) + list(_IV[:4]) + [counter & 0xFFFFFFFF, (counter >> 32) & 0xFFFFFFFF, block_len, flags]
    m = list(words)
    for _ in range(7):
        _g(st, 0, 4, 8, 12, m[0], m[1])
        _g(st, 1, 5, 9, 13, m[2], m[3])
        _g(st, 2, 6, 10, 14, m[4], m[5])
        _g(st, 3, 7, 11, 15, m[6], m[7])
        _g(st, 0, 5, 10, 15, m[8], m[9])
        _g(st, 1, 6, 11, 12, m[10], m[11])
        _g(st, 2, 7, 8, 13, m[12], m[13])
        _g(st, 3, 4, 9, 14, m[14], m[15])
        m = [m[i] for i in _PERMUTATION]
    return [st[i] ^ st[i + 8] for i in range(8)] + [st[i + 8] ^ cv[i] for i in range(8)]


def _words(block: bytes) -> list[int]:
    block = block + bytes(64 - len(block))
    return [int.from_bytes(block[i:i + 4], "little") for i in range(0, 64, 4)]


def _chunk_output(chunk: bytes, counter: int):
    blocks = [chunk[i:i + 64] for i in range(0, len(chunk), 64)] or [b""]
    cv = list(_IV)
    for i, block in enumerate(blocks[:-1]):
        cv = _compress(cv, _words(block), counter, 64, _CHUNK_START if i == 0 else 0)[:8]
    flags = _CHUNK_END | (_CHUNK_START if len(blocks) == 1 else 0)
    return cv, _words(blocks[-1]), counter, len(blocks[-1]), flags


def blake3(data: bytes, length: int = 32) -> bytes:
    """BLAKE3's extended output: its first 32 bytes are the digest."""
    chunks = [data[i:i + 1024] for i in range(0, len(data), 1024)] or [b""]
    stack: list[list[int]] = []
    for index, chunk in enumerate(chunks[:-1]):
        cv = _compress(*_chunk_output(chunk, index))[:8]
        total = index + 1
        while total & 1 == 0:
            cv = _compress(list(_IV), stack.pop() + cv, 0, 64, _PARENT)[:8]
            total >>= 1
        stack.append(cv)
    out = _chunk_output(chunks[-1], len(chunks) - 1)
    while stack:
        left = stack.pop()
        right = _compress(*out)[:8]
        out = (list(_IV), left + right, 0, 64, _PARENT)
    stream = b""
    block = 0
    while len(stream) < length:
        words = _compress(out[0], out[1], block, out[3], out[4] | _ROOT)
        stream += b"".join(w.to_bytes(4, "little") for w in words)
        block += 1
    return stream[:length]


# -- design §2's hashed target ------------------------------------------------------------------

DOMAIN = b"ic-workflow-public-target-v1\0"


def hashed_target(curve, a: int, h: int, r: int, seed: int):
    """v1's `public_hash_target` (src/bin/ic/workflow.rs), read past one
    word as design §2 says: x from the first W little-endian words of the
    extended output, masked to n bits; byte 8W chooses between the two
    lifts in ascending order of y; the cofactor clears it."""
    n = curve.n
    w = (n + 63) // 64
    for counter in range(1_000_000):
        data = (DOMAIN + n.to_bytes(4, "little") + bytes([a]) + (1).to_bytes(4, "little")
                + (1).to_bytes(8, "little") + seed.to_bytes(8, "little") + counter.to_bytes(8, "little"))
        stream = blake3(data, 8 * w + 1)
        x = int.from_bytes(stream[:8 * w], "little") & ((1 << n) - 1)
        lifts = sorted(curve.lift_all(x))
        if not lifts:
            continue
        raw = lifts[(stream[8 * w] & 1) % len(lifts)]
        target = curve.mul(h, raw)
        if target is None:
            continue
        assert curve.mul(r, target) is None
        return target, counter
    raise SystemExit("no hashed target within the attempt cap")


class Curve(v2.Curve):
    def lift_all(self, x: int) -> list[tuple[int, int]]:
        """Both points with abscissa x (one when x = 0), as (x, y)."""
        p = self.lift(x)
        if p is None:
            return []
        q = (p[0], p[1] ^ p[0])
        return [p] if q == p else [p, q]


# -- the instances' free action -------------------------------------------------------------------

def sqrt_mod(c: int, p: int) -> int | None:
    """A square root of c modulo the odd prime p (Tonelli–Shanks), or None."""
    c %= p
    if c == 0:
        return 0
    if pow(c, (p - 1) // 2, p) != 1:
        return None
    q, s = p - 1, 0
    while q % 2 == 0:
        q, s = q // 2, s + 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, cc, t, root = s, pow(z, q, p), pow(c, q, p), pow(c, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2, i = t2 * t2 % p, i + 1
        b = pow(cc, 1 << (m - i - 1), p)
        m, cc, t, root = i, b * b % p, t * b * b % p, root * b % p
    return root


def acts_freely(a: int, n: int, r: int) -> bool:
    """The eigenvalue of Frobenius on the subgroup of order r, a root of
    x^2 - mu x + 2 with lambda^n = 1, has order exactly n. For prime n that
    is lambda != 1."""
    assert all(n % d for d in range(2, int(n ** 0.5) + 1)), "the instances' degrees are prime"
    mu = 1 if a == 1 else -1
    s = sqrt_mod(mu * mu - 8, r)
    if s is None:
        return False
    inv2 = pow(2, -1, r)
    for root in ((mu + s) * inv2 % r, (mu - s) * inv2 % r):
        if pow(root, n, r) == 1:
            return root != 1
    return False


# -- the cases ---------------------------------------------------------------------------------------

def case(cid: str, purpose: str, files: dict, argv: list[str], expect: dict, timeout: int = TIMEOUT) -> dict:
    return {"id": cid, "step": STEP, "checks": purpose, "files": files, "argv": argv,
            "env": {"RAYON_NUM_THREADS": "1"}, "expect": expect, "timeout_s": timeout}


def price(params: str = "{tmp}/params.json", report: str = "{tmp}/report.json") -> list[str]:
    return ["price", "--params", params, "--json", "--out", report, "--repeats", "1", "--repeats-fast", "1"]


def check(params: str = "{tmp}/params.json", report: str = "{tmp}/report.json") -> list[str]:
    return ["check", "--params", params, "--json", "--out", report]


def known_log(slug: str, r: int) -> int:
    digest = hashlib.sha256(f"{LABEL}/{slug}/known_log".encode()).digest()
    return int.from_bytes(digest, "big") % (r - 1) + 1


def build() -> tuple[dict[str, dict | str], list[dict]]:
    files: dict[str, dict | str] = {}
    cases: list[dict] = []
    both = [{"pipeline": "kic", "admitted": True}, {"pipeline": "rho-koblitz", "admitted": True}]

    # C088-C093: kic at F0 on three-word fields.
    for cid, a, n, r, purpose in INSTANCES:
        f = v2.curve_id.find_irreducible_sparse(n)
        assert v2.irreducible(f) and 127 <= n <= 190, n
        curve = Curve(n, f, a, 1)
        order = v2.koblitz_order(a, n)
        h = order // r
        assert order == r * h and v2.prime_exact(r) and h % r, (a, n)
        assert acts_freely(a, n, r), (a, n, r)
        ids = v2.curve_id.binary_id(n, f, a, 1, order, end="-7")
        ident, slug = ids["icv1"], ids["slug"]
        gen = curve.subgroup_point(f"{LABEL}/{ident}/generator", h, r)
        k = known_log(ident, r)
        name = f"{cid}: {slug}, {purpose}"
        fname = f"{cid}-kic-three-word-n{n}.json"
        files[fname] = v2.document(name, curve, {"form": "koblitz", "a": a}, r, h, gen,
                                   {"known_log": str(k)}, "auto", RHO_SEED, routed=True)
        # Measurement 5's public target: a point hashed from a public
        # label, whose logarithm nobody knows.
        t001 = curve.subgroup_point(f"{LABEL}/{ident}/T001", h, r)
        files[f"{cid}-T001-n{n}.json"] = v2.document(
            f"{cid}-T001: {slug}, a public target for B4's measurement 5", curve, {"form": "koblitz", "a": a},
            r, h, gen, {"point": v2.point_doc(t001)}, "auto", RHO_SEED, routed=True)
        cases.append(case(
            f"{cid}-kic-three-word-n{n}", f"B4's kic at F0 on a three-word field: {purpose}",
            {"params.json": {"copy": f"{{here}}/{fname}"}}, price(),
            {"exit": 0, "json_file": "{tmp}/report.json",
             "json_paths": {"status": "complete", "fidelity": "F0", "result.scalar": str(k),
                            "result.verified": True, "result.known_answer": True, "all_verified": True,
                            "ic_and_rho_agree": True, "curve_id.icv1": ident, "ic.words": 3},
             "json_contains": {"route.considered": both}}))

    # C094, C095: sameness end to end. `--kic-multi` runs the multi-word
    # pipeline, at three words, where a narrower one would run.
    compared = ["counts", "certificates.ic.scalar", "certificates.rho.scalar", "rho_counts", "all_verified",
                "ic_and_rho_agree"]
    cases.append(case(
        "C094-smoke-multi-word-is-one-word", "the multi-word pipeline on the smoke row's translation gives "
        "the one-word pipeline's counts, logarithm and certificates",
        {"params.json": {"copy": SMOKE}}, price() + ["--kic-multi"],
        {"exit": 0, "json_file": "{tmp}/report.json",
         "json_paths": {"status": "complete", "ic.words": 3, "all_verified": True},
         "same_outputs_as": {"argv": price("{tmp}/params.json", "{tmp}/narrow.json"),
                             "json_file": "{tmp}/narrow.json", "paths": compared}},
        timeout=300))
    files[C080_DOC] = (B3B / "params" / C080_DOC).read_text()
    cases.append(case(
        "C095-n79-multi-word-is-two-word", "C094 on B3b's C080 document (n = 79): the multi-word pipeline gives "
        "the two-word pipeline's counts, logarithm and certificates",
        {"params.json": {"copy": f"{{here}}/{C080_DOC}"}}, price() + ["--kic-multi"],
        {"exit": 0, "json_file": "{tmp}/report.json",
         "json_paths": {"status": "complete", "ic.words": 3, "all_verified": True},
         "same_outputs_as": {"argv": price("{tmp}/params.json", "{tmp}/narrow.json"),
                             "json_file": "{tmp}/narrow.json", "paths": compared}},
        timeout=900))

    # C096: C055's and C086's successor (both `until: B4`). n = 131 now
    # passes the field gate; r = 2^129 meets the scalar gate, which B4
    # keeps (design §2) and B7b lifts.
    scalar_gate = [{"pipeline": p, "admitted": False, "gate": "scalar-wider-than-127-bits"}
                   for p in ("kic", "rho-koblitz")]
    cases.append(case(
        "C096-challenge-past-the-field-gate", "C055's and C086's successor from B4: the challenge file is valid; "
        "n = 131 passes the field gate, and kic and rho-koblitz refuse r = 2^129 as wider than 127 bits",
        {"params.json": {"copy": CHALLENGE}}, check(),
        {"exit": 0, "json_file": "{tmp}/report.json", "json_paths": {"status": "checks_passed"},
         "json_contains": {"checks": [{"code": "order-composite", "status": "pass", "exact": False}],
                           "disclosures": [{"code": "primality-screen"}],
                           "route.considered": scalar_gate}},
        timeout=120))

    # C097: C056's successor (`until: B4`): v1's hashed target at n = 83.
    gate = json.loads((V2 / "params" / "gate-m83-T001.json").read_text())
    assert int(gate["field"]["modulus"], 16) == GATE_MODULUS and gate["field"]["degree"] == 83
    g_curve = Curve(83, GATE_MODULUS, gate["curve"]["a"], 1)
    g_r, g_h = int(gate["subgroup"]["order"]), int(gate["subgroup"]["cofactor"])
    target, counter = hashed_target(g_curve, gate["curve"]["a"], g_h, g_r, 1)
    cases.append(case(
        "C097-gate-rho-hashed-target", "C056's successor from B4: v1's hashed target past one word is derived "
        "by design §2's rule, equal to this generator's own derivation, and rho alone stops at its step cap",
        {"params.json": {"copy": GATE, "set": {"target": {"public_hash_seed": 1}, "method": {
            "solve": "rho", "fidelity": "F0",
            "rho": {"pipeline": "auto", "seed": RHO_SEED, "max_iterations": 1000}}}}},
        ["price", "--params", "{tmp}/params.json", "--json", "--out", "{tmp}/report.json"],
        {"exit": 1, "json_file": "{tmp}/report.json",
         "json_paths": {"status": "not_recovered", "result.verified": False, "rho.step_cap": 1000,
                        "rho.pipeline": "rho-koblitz",
                        "target.record.kind": "public_hash_to_curve_cofactor",
                        "target.record.public_hash_seed": 1, "target.record.public_hash_counter": counter,
                        "target.record.x": str(target[0]), "target.record.y": str(target[1])}},
        timeout=120))

    # C098: past nine words, the schema's degree range refuses first.
    f = v2.curve_id.find_irreducible_sparse(WIDE_N)
    assert v2.irreducible(f)
    curve = Curve(WIDE_N, f, WIDE_A, 1)
    order = v2.koblitz_order(WIDE_A, WIDE_N)
    h = order // WIDE_R
    assert order == WIDE_R * h and v2.prime_exact(WIDE_R) and h % WIDE_R
    gen = curve.subgroup_point(f"{LABEL}/n{WIDE_N}/generator", h, WIDE_R)
    slug = v2.curve_id.binary_id(WIDE_N, f, WIDE_A, 1, order, end="-7")["slug"]
    k = known_log(slug, WIDE_R)
    # The slug runs to 113 characters at this degree, so the name, which
    # schema v2 caps at 120, describes the curve instead.
    name = f"C098: the Koblitz curve with a = {WIDE_A} at n = {WIDE_N}, past nine words' limit of 574"
    assert len(name) <= 120 and len(slug) > 100
    files["C098-past-nine-words-n577.json"] = v2.document(
        name, curve, {"form": "koblitz", "a": WIDE_A}, WIDE_R, h, gen,
        {"known_log": str(k)}, "auto", RHO_SEED, routed=True)
    cases.append(case(
        "C098-past-nine-words", "past nine words: schema v2's binary degrees stop at 571, inside nine "
        "words' 574, so n = 577 is refused as degree-range before any gate",
        {"params.json": {"copy": "{here}/C098-past-nine-words-n577.json"}}, check(),
        {"exit": 2, "json_file": "{tmp}/report.json",
         "json_paths": {"status": "refused", "refusal.code": "degree-range", "refusal.class": "invalid"}},
        timeout=120))

    # C099: sect163k1, paired at F0: admitted by width, refused by budget.
    s = SECT163K1
    curve = Curve(163, s["modulus"], 1, 1)
    order = v2.koblitz_order(1, 163)
    assert order == s["r"] * s["h"] and curve.on_curve(s["g"]) and curve.mul(s["r"], s["g"]) is None
    files["C099-sect163k1-paired.json"] = v2.document(
        "C099: sect163k1 (SEC 2), paired at F0", curve, {"form": "koblitz", "a": 1}, s["r"], s["h"],
        s["g"], {"known_log": "1"}, "auto", RHO_SEED, routed=True)
    cases.append(case(
        "C099-sect163k1-past-the-field-gate", "sect163k1: valid, n = 163 passes the field gate at three words, "
        "and kic and rho-koblitz refuse r = 2^162 as wider than 127 bits",
        {"params.json": {"copy": "{here}/C099-sect163k1-paired.json"}}, check(),
        {"exit": 0, "json_file": "{tmp}/report.json", "json_paths": {"status": "checks_passed"},
         "json_contains": {"route.considered": scalar_gate}},
        timeout=120))

    # C102: a three-word subgroup both gates admit, too large for F0.
    # #K_0(GF(2^127)) = 4 * 1118452171 * WIDE_R95, the last a 95-bit
    # prime (its primality here is a screen: it is past the exact range
    # of v2's test).
    f = v2.curve_id.find_irreducible_sparse(127)
    curve = Curve(127, f, 0, 1)
    order = v2.koblitz_order(0, 127)
    h = order // R95
    assert order == R95 * h and h == 4 * 1118452171 and v2.miller_rabin(R95, PRIMES_40)
    assert (1 << 63) < R95 < (1 << 127)
    gen = curve.subgroup_point(f"{LABEL}/n127-r95/generator", h, R95)
    files["C102-n127-95-bit-subgroup.json"] = v2.document(
        f"C102: {v2.curve_id.binary_id(127, f, 0, 1, order, end='-7')['slug']}, its 95-bit subgroup, paired at F0",
        curve, {"form": "koblitz", "a": 0}, R95, h, gen,
        {"known_log": "12345"}, "auto", RHO_SEED, routed=True)
    cases.append(case(
        "C102-n127-over-budget", "a 95-bit subgroup at n = 127: kic and rho-koblitz admit it by field and scalar "
        "width, and F0's estimate exceeds the budget, so the paired price is refused as over budget",
        {"params.json": {"copy": "{here}/C102-n127-95-bit-subgroup.json"}}, price(),
        {"exit": 4, "json_file": "{tmp}/report.json",
         "json_paths": {"refusal.code": "over-budget", "refusal.class": "over_budget"},
         "json_contains": {"route.considered": both}},
        timeout=300))

    # C100, C101: v1's other rules past one word (design §2), at B3b's
    # smallest two-word instance, end to end. Their sameness with v1's
    # rules at one word is a unit test.
    files[C078_DOC] = (B3B / "params" / C078_DOC).read_text()
    for cid, purpose, edit in (
            ("C100-random-target-past-one-word", "v1's random target past one word: drawn by design §2's rule at "
             "n = 67, recovered and verified by both arms", {"target": {"random_seed": 7}}),
            ("C101-generator-rule-past-one-word", "v1's generator rule past one word: found by design §2's rule at "
             "n = 67, and the known logarithm recovered and verified by both arms",
             {"subgroup.generator": {"rule": "koblitz_search_v1"}})):
        cases.append(case(
            cid, purpose, {"params.json": {"copy": f"{{here}}/{C078_DOC}", "set": edit}}, price(),
            {"exit": 0, "json_file": "{tmp}/report.json",
             "json_paths": {"status": "complete", "fidelity": "F0", "result.verified": True,
                            "result.known_answer": True, "all_verified": True, "ic_and_rho_agree": True,
                            "ic.words": 2},
             "json_contains": {"route.considered": both}}))
    cases.sort(key=lambda c: c["id"])
    return files, cases


def selftest() -> None:
    """BLAKE3's published test vectors: the empty input, whose 131-byte
    extended output spans three output blocks; "abc"; and the vectors'
    1,025-byte input (bytes i mod 251), which spans two chunks."""
    assert blake3(b"abc").hex() == "6437b3ac38465133ffb63b75273a8db548c558465d79db03fd359c6cd5bd9d85"
    assert blake3(bytes(i % 251 for i in range(1025))).hex() == (
        "d00278ae47eb27b34faecf67b4fe263f82d5412916c1ffd97c8cb7fb814b8444")
    assert blake3(b"").hex() == "af1349b9f5f9a1a6a0404dea36dcc9499bcb25c9adc112b7cc9a93cae41f3262"
    assert blake3(b"", 131).hex() == (
        "af1349b9f5f9a1a6a0404dea36dcc9499bcb25c9adc112b7cc9a93cae41f3262e00f03e7b69af26b7faaf09fcd333050"
        "338ddfe085b8cc869ca98b206c08243a26f5487789e8f660afe6c99ef9e0c52b92e7393024a80459cf91f476f9ffdbda"
        "7001c22e159b402631f277ca96f2defdf1078282314e763699a31c5363165421cce14d")


def texts() -> dict[str, str]:
    selftest()
    files, cases = build()
    assert all(len(doc["name"]) <= 120 for doc in files.values() if isinstance(doc, dict)), "schema v2: name"
    out = {f"params/{name}": (doc if isinstance(doc, str) else json.dumps(doc, indent=1) + "\n")
           for name, doc in files.items()}
    out["cases.json"] = json.dumps({
        "suite": "ic tool programme conformance suite v2, B4's cases",
        "design": "research/ic_tool_program/design/multi-word.md; "
                  "research/ic_tool_program/rounds/B4-multi-word/PROTOCOL.md",
        "includes": "../v1/cases.json (step B0), ../v2/cases.json (B1) and the later steps' sets, run first, "
                    "with the until rule",
        "rules": [
            "The rules of ../v2/cases.json hold. {here} is this directory's params/; {cases} is still "
            "../v2/params/. ../run.py runs these after the earlier steps' cases.",
            "The until rule: C031, C055, C056 and C086 name B4 as their `until`, so from B4 they are not run. "
            "C096 succeeds C055 and C086, and C097 succeeds C056.",
            "`--kic-multi` is B4's test switch: it runs the multi-word pipeline, at three words, where the one- "
            "or two-word pipeline would run. C094 and C095 compare them through `same_outputs_as`.",
            "params/C078-kic-two-word-n67.json and params/C080-kic-two-word-n79.json are B3b's files, copied byte "
            "for byte for C100, C101 and C095.",
            "params/C088-T001-n127.json to params/C093-T001-n179.json are no case's: they are B4's measurement "
            "5's public targets, frozen here with the cases' documents.",
        ],
        "generator": {"path": "research/ic_tool_program/conformance/v2-b4/make_cases.py",
                      "sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        "cases": cases,
    }, indent=1, ensure_ascii=False) + "\n"
    out["SHA256SUMS"] = "".join(f"{hashlib.sha256(t.encode()).hexdigest()}  {rel}\n"
                                for rel, t in sorted(out.items()))
    return out


def main() -> None:
    out = texts()
    if "--check" in sys.argv:
        bad = [rel for rel, t in out.items() if not (HERE / rel).exists() or (HERE / rel).read_text() != t]
        print(json.dumps({"files": len(out), "mismatches": bad}, indent=1))
        raise SystemExit(1 if bad else 0)
    if (HERE / "cases.json").exists():
        raise SystemExit("cases.json exists; B4's cases are frozen (use --check)")
    for rel, t in out.items():
        path = HERE / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "x") as fh:
            fh.write(t)


if __name__ == "__main__":
    main()
