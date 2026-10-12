"""
ecmin.py -- toy elliptic-curve arithmetic and the Abdullah-Mahalanobis-Mallick
"minors" relation matrix, in pure Python (no numpy/sympy required).

Background
----------
For an elliptic curve E/F_p in short Weierstrass form y^2 = x^3 + a x + b, the
Riemann-Roch space L(k.O) (functions with a pole of order <= k only at the point
at infinity O) has dimension k and a monomial basis ordered by pole order:

    pole order:  0   2   3   4    5    6    7    8     9   ...
    function  :  1   x   y   x^2  xy   x^3  x^2 y x^4   x^3 y ...

(pole order of x is 2, of y is 3; y^2 is reduced via the curve equation, so each
function is x^i y^j with j in {0,1}, sorted by 2i+3j.)

Key fact (the linear-algebra form of Semaev's relation):
    k points P_1,...,P_k on E sum to O in the group
      <=>  the k x k matrix  [ phi_j(P_i) ]_{i,j=1..k}   is SINGULAR over F_p,
where phi_1,...,phi_k are the first k basis functions of L(k.O).

The Abdullah-Mahalanobis-Mallick method builds the N x N matrix
    M[i][j] = phi_j(R_i),   R_i = a_i P + b_i Q,
and searches for a vanishing minor. A vanishing k-subset gives
    sum_{i in S} (a_i + m b_i) = 0  (mod n)   =>   m = -(sum a_i)/(sum b_i),
solving the ECDLP. The leading principal k x k minors are the "initial minors";
they are exactly the pivots of an LU / Gaussian elimination of M (Schur
complements), so a zero pivot at step k certifies that the first k rows sum to O.
"""

import random


# ----------------------------------------------------------------------------
# Field / curve arithmetic
# ----------------------------------------------------------------------------

def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    ls = pow(a, (p - 1) // 2, p)
    return -1 if ls == p - 1 else 1


def curve_order(a, b, p):
    """#E(F_p) by summing Legendre symbols. O(p); only for toy p."""
    total = p + 1  # account for O plus the contribution sum
    for x in range(p):
        total += legendre((x * x * x + a * x + b) % p, p)
    return total


class Curve:
    """y^2 = x^3 + a x + b over F_p.  O represented as None."""

    def __init__(self, a, b, p):
        self.a = a % p
        self.b = b % p
        self.p = p

    def is_on(self, P):
        if P is None:
            return True
        x, y = P
        return (y * y - (x * x * x + self.a * x + self.b)) % self.p == 0

    def add(self, P, Q):
        p = self.p
        if P is None:
            return Q
        if Q is None:
            return P
        x1, y1 = P
        x2, y2 = Q
        if x1 == x2 and (y1 + y2) % p == 0:
            return None
        if P == Q:
            num = (3 * x1 * x1 + self.a) % p
            den = (2 * y1) % p
        else:
            num = (y2 - y1) % p
            den = (x2 - x1) % p
        lam = (num * pow(den, p - 2, p)) % p
        x3 = (lam * lam - x1 - x2) % p
        y3 = (lam * (x1 - x3) - y1) % p
        return (x3, y3)

    def mul(self, k, P):
        R = None
        k %= (self.p * 4)  # generous; callers pass k < group order anyway
        Q = P
        while k > 0:
            if k & 1:
                R = self.add(R, Q)
            Q = self.add(Q, Q)
            k >>= 1
        return R

    def find_generator_of_prime_order(self, n):
        """Return a point of order exactly n (n must divide #E)."""
        while True:
            x = random.randrange(self.p)
            rhs = (x * x * x + self.a * x + self.b) % self.p
            if legendre(rhs, self.p) < 0:
                continue
            y = pow(rhs, (self.p + 1) // 4, self.p) if self.p % 4 == 3 else _tonelli(rhs, self.p)
            if y is None:
                continue
            P = (x, y)
            # cofactor-clear handled by caller; here assume #E = n (prime order)
            if self.mul(n, P) is None and P is not None:
                return P


def _tonelli(a, p):
    a %= p
    if a == 0:
        return 0
    if legendre(a, p) != 1:
        return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    # Tonelli-Shanks
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while legendre(z, p) != -1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        t2, i = t, 0
        while t2 != 1:
            t2 = (t2 * t2) % p
            i += 1
        bexp = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, (bexp * bexp) % p, (t * bexp * bexp) % p, (r * bexp) % p
    return r


# ----------------------------------------------------------------------------
# Riemann-Roch monomial basis and the relation matrix
# ----------------------------------------------------------------------------

def rr_basis_exponents(N):
    """First N functions of L(k.O), as (i, j) meaning x^i y^j, sorted by pole
    order 2i + 3j.  j in {0,1}."""
    cands = []
    bound = 2 * N + 4
    for i in range(bound):
        for j in (0, 1):
            order = 2 * i + 3 * j
            cands.append((order, i, j))
    cands.sort()
    return [(i, j) for (_, i, j) in cands[:N]]


def eval_phi(point, i, j, p):
    """phi(P) = x^i y^j.  At O the affine functions blow up; callers never pass
    O as a row (we only use affine points)."""
    x, y = point
    return (pow(x, i, p) * pow(y, j, p)) % p


def build_matrix(points, p, ncols=None):
    """Rows = affine points, columns = RR basis functions phi_1..phi_ncols."""
    n = len(points)
    ncols = ncols or n
    basis = rr_basis_exponents(ncols)
    M = [[eval_phi(P, i, j, p) for (i, j) in basis] for P in points]
    return M


# ----------------------------------------------------------------------------
# Determinant / rank over F_p
# ----------------------------------------------------------------------------

def submatrix_singular(M, rows, cols, p):
    A = [[M[r][c] for c in cols] for r in rows]
    return det_modp(A, p) == 0


def det_modp(A, p):
    A = [row[:] for row in A]
    n = len(A)
    det = 1
    for col in range(n):
        piv = None
        for r in range(col, n):
            if A[r][col] % p != 0:
                piv = r
                break
        if piv is None:
            return 0
        if piv != col:
            A[col], A[piv] = A[piv], A[col]
            det = (-det) % p
        inv = pow(A[col][col], p - 2, p)
        det = (det * A[col][col]) % p
        for r in range(col + 1, n):
            f = (A[r][col] * inv) % p
            if f:
                for c in range(col, n):
                    A[r][c] = (A[r][c] - f * A[col][c]) % p
    return det % p


def leading_minor_zero_pivots(M, p):
    """Gaussian elimination WITHOUT row swaps on the leading square block.
    Returns the index k (1-based) of the first step whose leading principal
    minor vanishes, i.e. the first prefix of rows that is dependent in L(kO);
    or None if all leading principal minors up to size = min(rows,cols) are
    nonzero.

    NOTE: we deliberately do not pivot, because the 'initial minor' / Schur
    recursion is about the *leading principal* minors of the rows in their given
    (random) order.  det(leading m x m) = prod of first m diagonal pivots.
    """
    rows = len(M)
    cols = len(M[0])
    size = min(rows, cols)
    A = [row[:] for row in M]
    for k in range(size):
        if A[k][k] % p == 0:
            # leading (k+1) principal minor is zero (given previous nonzero)
            return k + 1
        inv = pow(A[k][k], p - 2, p)
        for r in range(k + 1, rows):
            f = (A[r][k] * inv) % p
            if f:
                for c in range(k, cols):
                    A[r][c] = (A[r][c] - f * A[k][c]) % p
    return None
