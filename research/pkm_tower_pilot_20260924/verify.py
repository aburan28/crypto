#!/usr/bin/env python3
"""Independent cross-check of every pilot row against exhaustive search.

`AGENTS.md` section 6 does not accept a new oracle that has not been
cross-checked against the one it replaces on every input of a run. For each
row of the pilot this recomputes, without F4, how many tuples of `V^m` solve
the row's summation system:

- m = 2: ordered pairs (x1, x2) in V x V with S3(x1, x2, xR) = 0;
- m = 3: ordered triples with S4(x1, x2, x3, xR) = 0, where
  S4 = Res_U(S3(x1, x2, U), S3(x3, xR, U)). This is the chain system with its
  free unknown eliminated, so it counts solutions whose partial sum is only
  defined over F_{p^2}, exactly as the F_p ideal does.
- m = 4: ordered quadruples with S5 = Res_U(S4(x1, x2, x3, U), S3(U, x4, xR)) = 0,
  the chain with both free unknowns eliminated, as a Sylvester determinant.

It then checks the only thing that can be checked without trusting F4: a
refuted system (`inconsistent`) has no solution on the grid, and a system F4
did not refute has at least one. The tower's `V` is taken from the row: `V` for
isogeny rows, rebuilt from `g, zeta` (Kummer) or `u, zeta` (Dickson) otherwise.
Rows of the null, naive and ladder controls are skipped, because their
systems are not the summation system.

    python3 verify.py runs/*.jsonl
"""

import json
import sys


def s3(x1, x2, x3, a, b, p):
    d = (x1 - x2) % p
    s = (x1 + x2) % p
    pr = x1 * x2 % p
    t1 = d * d % p * x3 % p * x3 % p
    inner = (s * ((pr + a) % p) + 2 * b) % p
    t2 = (p - 2) * inner % p * x3 % p
    q = (pr - a) % p
    t3 = (q * q - 4 * b * s) % p
    return (t1 + t2 + t3) % p


def s3_coeffs(x1, x2, a, b, p):
    """S3(x1, x2, U) as A U^2 + B U + C."""
    d = (x1 - x2) % p
    s = (x1 + x2) % p
    pr = x1 * x2 % p
    A = d * d % p
    B = (p - 2) * ((s * ((pr + a) % p) + 2 * b) % p) % p
    q = (pr - a) % p
    C = (q * q - 4 * b * s) % p
    return A, B, C


def res2(f, g, p):
    """Resultant of two quadratics A U^2 + B U + C."""
    A, B, C = f
    A2, B2, C2 = g
    u = (A * C2 - A2 * C) % p
    v = (A * B2 - A2 * B) % p
    w = (B * C2 - B2 * C) % p
    return (u * u - v * w) % p


def pmul(f, g, p):
    """Product of two polynomials in U, as coefficient lists (index = power)."""
    out = [0] * (len(f) + len(g) - 1)
    for i, a in enumerate(f):
        for j, b in enumerate(g):
            out[i + j] = (out[i + j] + a * b) % p
    return out


def psub(f, g, p):
    n = max(len(f), len(g))
    f = f + [0] * (n - len(f))
    g = g + [0] * (n - len(g))
    return [(a - b) % p for a, b in zip(f, g)]


def s4_in_u(x1, x2, x3, a, b, p):
    """S4(x1, x2, x3, U) = Res_w(S3(x1, x2, w), S3(w, x3, U)), as a polynomial in U.

    S3(x1, x2, w) = A w^2 + B w + C with constant coefficients; S3(w, x3, U),
    read as a quadratic in w, has coefficients that are quadratics in U. The
    two-quadratic resultant formula applies with polynomial coefficients.
    """
    A, B, C = s3_coeffs(x1, x2, a, b, p)
    # S3(x3, U, w) as a quadratic in w: coefficients as polynomials in U (low to high).
    A2 = [x3 * x3 % p, (p - 2) * x3 % p, 1]
    B2 = [(p - 2) * ((a * x3 + 2 * b) % p) % p, (p - 2) * ((x3 * x3 + a) % p) % p, (p - 2) * x3 % p]
    C2 = [(a * a - 4 * b * x3) % p, (p - 2) * ((a * x3 + 2 * b) % p) % p, x3 * x3 % p]
    u = psub([A * c % p for c in C2], [C * c % p for c in A2], p)
    v = psub([A * c % p for c in B2], [B * c % p for c in A2], p)
    w = psub([B * c % p for c in C2], [C * c % p for c in B2], p)
    return psub(pmul(u, u, p), pmul(v, w, p), p)


def res_zero(f, g, p):
    """Whether Res(f, g) = 0 for polynomials in U (coefficient lists, low to high),
    as the Sylvester determinant at their formal degrees."""
    df, dg = len(f) - 1, len(g) - 1
    n = df + dg
    rows = []
    for i in range(dg):
        rows.append([0] * i + f[::-1] + [0] * (n - df - 1 - i))
    for i in range(df):
        rows.append([0] * i + g[::-1] + [0] * (n - dg - 1 - i))
    # Gaussian elimination mod p; the determinant is zero iff some column has no pivot.
    for c in range(n):
        piv = next((r for r in range(c, n) if rows[r][c] % p), None)
        if piv is None:
            return True
        rows[c], rows[piv] = rows[piv], rows[c]
        inv = pow(rows[c][c], p - 2, p)
        for r in range(c + 1, n):
            if rows[r][c] % p:
                k = rows[r][c] * inv % p
                rows[r] = [(x - k * y) % p for x, y in zip(rows[r], rows[c])]
    return False


def tower_v(row):
    info = row["tower"]
    p = row["p"]
    t = row["t"]
    if "v" in info:
        return info["v"]
    if row["kind"] == "kummer":
        g, z = info["g"], info["zeta"]
        return [g * pow(z, k, p) % p for k in range(1 << t)]
    if row["kind"] == "dickson":
        u, z = info["u"], info["zeta"]
        out = []
        for k in range(1 << t):
            w = u * pow(z, k, p) % p
            out.append((w + pow(w, p - 2, p)) % p)
        return out
    return None


def count(row, v):
    p = row["p"]
    a, b = row["curve"]["a"], row["curve"]["b"]
    xr = row["x_r"]
    if row["m"] == 2:
        return sum(1 for x1 in v for x2 in v if s3(x1, x2, xr, a, b, p) == 0)
    if row["m"] == 3:
        tails = {x3: s3_coeffs(x3, xr, a, b, p) for x3 in v}
        n = 0
        for x1 in v:
            for x2 in v:
                f = s3_coeffs(x1, x2, a, b, p)
                for x3 in v:
                    if res2(f, tails[x3], p) == 0:
                        n += 1
        return n
    if row["m"] == 4:
        # S5 = Res_U(S4(x1, x2, x3, U), S3(U, x4, x_R)) and S3(U, x4, x_R) = S3(x4, x_R, U).
        tails = {}
        for x4 in v:
            A, B, C = s3_coeffs(x4, xr, a, b, p)
            tails[x4] = [C, B, A]
        n = 0
        for x1 in v:
            for x2 in v:
                for x3 in v:
                    s4 = s4_in_u(x1, x2, x3, a, b, p)
                    for x4 in v:
                        if res_zero(s4, tails[x4], p):
                            n += 1
        return n
    return None


def main(paths):
    checked = agreed = skipped = bounded = 0
    bad = []
    seen = set()
    for path in paths:
        with open(path) as fh:
            for line in fh:
                row = json.loads(line)
                if "summary" in row or row["control"] != "tower" or row["timed_out"]:
                    continue
                # A run that stopped at its degree bound with pairs left (a
                # confirmation run, note section 11.4) claims nothing about
                # solutions: it neither refuted nor finished.
                if row.get("pairs_above_bound", 0) > 0 and not row["inconsistent"]:
                    bounded += 1
                    continue
                # A system measured by two runs is checked once (see analyze.py).
                key = (row["p"], row["kind"], row["m"], row["t"], row["target"],
                       row["target_index"], row["x_r"], row["curve"]["a"], row["curve"]["b"],
                       json.dumps(row["tower"], sort_keys=True), row.get("engine", "f4_fp"))
                if key in seen:
                    continue
                seen.add(key)
                v = tower_v(row)
                if v is None:
                    skipped += 1
                    continue
                n = count(row, v)
                checked += 1
                ok = (n == 0) == row["inconsistent"]
                if row["target"] == "planted":
                    ok = ok and n >= 1
                if ok:
                    agreed += 1
                else:
                    bad.append((path, row.get("engine", "f4_fp"), row["kind"], row["m"], row["N"],
                                row["target"], n, row["inconsistent"]))
    print(f"checked {checked} tower rows against exhaustive search; {agreed} agree; {skipped} skipped (no V); "
          f"{bounded} stopped at their degree bound, with no verdict to check.")
    for b in bad:
        print("DISAGREE", b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
