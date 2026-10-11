#!/usr/bin/env python3
"""Turn gap_pairs.rs result files into ecbench specs (prime_explicit curves).

Every curve in a pair has prime order N, so any affine point generates.
The generator is the point with the smallest x having a square right-hand
side, and the smaller of its two y's: deterministic, so the spec is
re-derivable from the result file alone.
"""
import json, re, sys

def sqrt_mod(n, p):
    n %= p
    if n == 0: return 0
    if pow(n, (p - 1) // 2, p) != 1: return None
    if p % 4 == 3: return pow(n, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    m, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while t != 1:
        i, tt = 0, t
        while tt != 1: tt = tt * tt % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r

def generator(p, a, b):
    x = 0
    while True:
        y = sqrt_mod((x * x * x + a * x + b) % p, p)
        if y is not None and y != 0:
            return x, min(y, p - y)
        x += 1

def parse(path):
    txt = open(path).read()
    q = int(re.search(r"q = (\d+)", txt).group(1))
    ell = int(re.search(r"ell = (\d+)", txt).group(1))
    N = int(re.search(r"N = #E = (\d+)", txt).group(1))
    curves = []
    for m in re.finditer(r"(crater|floor) E(\d+): a = (\d+) b = (\d+)", txt):
        curves.append((m.group(1), int(m.group(3)), int(m.group(4))))
    return q, ell, N, curves

def spec(path, arms, rounds, timeout, label_suffix="", targets=4, target_seed=20261010):
    q, ell, N, curves = parse(path)
    cs = []
    for level, a, b in curves:
        gx, gy = generator(q, a, b)
        cs.append({"kind": "prime_explicit", "p": q, "a": a, "b": b,
                   "group_order": N, "r": N, "gx": gx, "gy": gy})
    return {
        "schema": "ecbench.spec/v1",
        "label": f"conductor gap ell={ell}: crater (f_E=1) vs floor (f_E={ell}) over GF(q), q {q.bit_length()} bits{label_suffix}",
        "description": (f"CM pairs with D=-19: t^2-4q = -19*ell^2, #E = N = {N} prime on every curve. "
                        f"Curve 1 is the crater (j = -884736, End = O_K); the rest are floor curves "
                        f"(End has conductor ell = {ell}). Same q, same N, same group structure; only the "
                        f"endomorphism conductor differs. Paired cost(floor)/cost(crater) is the measurement. "
                        f"Source: research/conductor_gap_cm_pairs_20261010/results/{path.split('/')[-1]}."),
        "workloads": {"curves": cs, "targets_per_curve": targets, "target_seed": target_seed},
        "arms": arms,
        "measurement": {"rounds": rounds, "warmup": 1, "order": "alternate", "seed": 1,
                        "isolation_required": "L2", "timeout_seconds": timeout},
    }

GENERIC = [
    {"name": "rho-neg", "role": "reference", "method": {"id": "rho.negation"}},
    {"name": "bsgs-neg", "role": "candidate", "method": {"id": "bsgs.negation"}},
    {"name": "kangaroo", "role": "candidate", "method": {"id": "kangaroo.vow"}},
    {"name": "rho-neg-aa", "role": "control", "method": {"id": "rho.negation"}},
]
IC = [{"name": "ic-mitm", "role": "candidate", "method": {"id": "ic.pipeline", "params": {
        "factor_base": "prime-abscissa:size=64", "oracle": "mitm:negation_folded=1"}}}]

if __name__ == "__main__":
    out = sys.argv[1]
    json.dump(spec("results/ell_263_q19.txt", GENERIC + IC, 5, 300), open(f"{out}/gap-ell263-q19.json", "w"), indent=1)
    json.dump(spec("results/ell_8209_q29.txt", GENERIC, 5, 600), open(f"{out}/gap-ell8209-q29.json", "w"), indent=1)
    json.dump(spec("results/ell_131101_q37.txt", GENERIC, 3, 1800), open(f"{out}/gap-ell131101-q37.json", "w"), indent=1)
    # Follow-up after the first three sessions (results/crater_vs_floor.md): the
    # index-calculus arm at ell = 263 showed floor/crater = 0.974 [0.953, 0.992]
    # on 4 crater workloads. PROTOCOL.md needs the same sign at a second ell
    # with the control arm at 1, so: 16 independent targets per curve, fresh
    # target seed, ic + reference + control only, at both ell.
    ICONLY = [{"name": "rho-neg", "role": "reference", "method": {"id": "rho.negation"}},
              {"name": "rho-neg-aa", "role": "control", "method": {"id": "rho.negation"}}] + IC
    json.dump(spec("results/ell_263_q19.txt", ICONLY, 3, 600, " (ic follow-up, 16 targets)", 16, 20261011),
              open(f"{out}/gap-ell263-q19-ic16.json", "w"), indent=1)
    json.dump(spec("results/ell_8209_q29.txt", ICONLY, 3, 1800, " (ic follow-up, 16 targets)", 16, 20261011),
              open(f"{out}/gap-ell8209-q29-ic16.json", "w"), indent=1)
    print("wrote 5 specs")
