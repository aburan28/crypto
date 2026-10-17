"""Independent integer-field replay of the frozen second-generation line screen."""

import hashlib
import json
from collections import Counter
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
plane = json.loads((HERE / "type_i_plane_quartic.json").read_text())
p = int(plane["p"])
assert p % 4 == 3
terms = [(tuple(map(int, key.split(","))), int(value) % p)
         for key, value in plane["quartic_coefficients"].items()]
assert all(sum(exponents) == 4 for exponents, _ in terms)


def points(name, skip_first=False):
    rows = [tuple(map(int, line.split())) for line in
            (HERE / name).read_text().splitlines()]
    return [row[1:] if skip_first else row for row in rows]


old = points("type_i_line_base_points.tsv", True)
new = points("type_i_line_new_points.tsv")
pairs = points("type_i_line_new_point_pairs.tsv")
assert len(old) == 141 and len(new) == len(pairs) == 9784
selected = old + new[:64]
existing = set(old + new)
assert len(existing) == 9925
original = set()
for k in range(32):
    i, j = pairs[2 * k]
    assert pairs[2 * k + 1] == (i, j)
    original.add(tuple(sorted((old[i - 1], old[j - 1],
                               new[2 * k], new[2 * k + 1]))))
assert len(original) == 32


def polynomial_for_line(m, b):
    f = [0] * 5
    for (ex, ev, _), c in terms:
        for j in range(ev + 1):
            f[ex + j] = (f[ex + j] + c * comb(ev, j)
                         * pow(m, j, p) * pow(b, ev - j, p)) % p
    return f


counts = dict(pairs=0, split=0, nonsplit=0, ramified=0, degree_other=0,
              one_existing=0, closed=0, original=0, novel=0)
novel_witnesses = []
original_multiplicity = Counter()
for i in range(141, 205):
    x1, y1 = selected[i]
    for j in range(i):
        x2, y2 = selected[j]
        assert x1 != x2
        counts["pairs"] += 1
        m = ((y2 - y1) * pow((x2 - x1) % p, -1, p)) % p
        b = (y1 - m * x1) % p
        f = polynomial_for_line(m, b)
        s = (x1 + x2) % p
        pr = x1 * x2 % p
        q2 = f[4]
        q1 = (f[3] + s * q2) % p
        q0 = (f[2] + s * q1 - pr * q2) % p
        assert (pr * q0 - f[0]) % p == 0
        assert (-s * q0 + pr * q1 - f[1]) % p == 0
        if q2 == 0:
            counts["degree_other"] += 1
            continue
        disc = (q1 * q1 - 4 * q2 * q0) % p
        if disc == 0:
            counts["ramified"] += 1
            continue
        if pow(disc, (p - 1) // 2, p) != 1:
            counts["nonsplit"] += 1
            continue
        counts["split"] += 1
        sd = pow(disc, (p + 1) // 4, p)
        assert sd * sd % p == disc
        inv = pow(2 * q2 % p, -1, p)
        roots = [(-q1 + sd) * inv % p, (-q1 - sd) * inv % p]
        assert roots[0] != roots[1]
        residual = [(x, (m * x + b) % p) for x in roots]
        hits = sum(Q in existing for Q in residual)
        if hits == 1:
            counts["one_existing"] += 1
        if hits == 2:
            counts["closed"] += 1
            quad = tuple(sorted((selected[i], selected[j], *residual)))
            if quad in original:
                counts["original"] += 1
                original_multiplicity[quad] += 1
            else:
                counts["novel"] += 1
                novel_witnesses.append([i + 1, j + 1, quad])
assert counts == dict(pairs=11040, split=5636, nonsplit=5404, ramified=0,
                      degree_other=0, one_existing=0, closed=160,
                      original=160, novel=0)
assert not novel_witnesses
assert len(original_multiplicity) == 32
assert set(original_multiplicity.values()) == {5}
gp_log = (HERE / "type_i_line_second_gen_closure.log").read_text()
assert "LINE_SECOND_GEN_PASS pairs=11040 split=5636 nonsplit=5404 ramified=0 degree_other=0 one_existing=0 closed=160 original=160 novel=0" in gp_log
assert "LINE_SECOND_GEN_NOVEL_WITNESSES=[]" in gp_log
out = {
    "status": "LINE_SECOND_GEN_INDEPENDENT_PASS",
    "counts": counts,
    "old_points": len(old),
    "selected_new_points": 64,
    "full_frozen_points": len(existing),
    "original_quartets": len(original),
    "original_representations_per_quartet": 5,
    "input_sha256": {name: sha(HERE / name) for name in
                     ("type_i_plane_quartic.json",
                      "type_i_line_base_points.tsv",
                      "type_i_line_new_points.tsv",
                      "type_i_line_new_point_pairs.tsv",
                      "type_i_line_second_gen_input.gp",
                      "line_second_gen_closure.gp",
                      "type_i_line_second_gen_closure.log")},
    "source_sha256": sha(__file__),
}
(HERE / "type_i_line_second_gen_verification.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print("LINE_SECOND_GEN_INDEPENDENT_PASS", json.dumps(counts, sort_keys=True))
