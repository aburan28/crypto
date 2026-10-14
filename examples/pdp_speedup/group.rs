use crypto_lib::cryptanalysis::semaev_decomp::Gf2;

/// Solver for the Artin–Schreier equation `z² + z = c` over the one-word
/// field, by Gaussian elimination on the `F₂`-linear map `z ↦ z² + z`.
#[derive(Clone, Debug)]
pub struct ArtinSchreier {
    /// Echelon rows `(value, preimage)` indexed by leading bit of `value`.
    rows: Vec<Option<(u64, u64)>>,
}

impl ArtinSchreier {
    pub fn new(gf: &Gf2) -> Self {
        let n = gf.n;
        let mut rows: Vec<Option<(u64, u64)>> = vec![None; n as usize];
        for i in 0..n {
            let e = 1u64 << i;
            let (mut v, mut w) = (gf.sqr(e) ^ e, e);
            while v != 0 {
                let b = 63 - v.leading_zeros() as usize;
                match rows[b] {
                    Some((pv, pw)) => {
                        v ^= pv;
                        w ^= pw;
                    }
                    None => {
                        rows[b] = Some((v, w));
                        break;
                    }
                }
            }
        }
        Self { rows }
    }

    /// A `z` with `z² + z = c`, or `None` when `Tr(c) = 1`.  The other
    /// solution is `z + 1`.
    pub fn solve(&self, c: u64) -> Option<u64> {
        let (mut v, mut w) = (c, 0u64);
        while v != 0 {
            let b = 63 - v.leading_zeros() as usize;
            let (pv, pw) = self.rows[b]?;
            v ^= pv;
            w ^= pw;
        }
        Some(w)
    }
}

/// A point of `y² + xy = x³ + x² + 1` over the one-word field.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Pt {
    Inf,
    Aff(u64, u64),
}

pub fn neg(p: Pt) -> Pt {
    match p {
        Pt::Inf => Pt::Inf,
        Pt::Aff(x, y) => Pt::Aff(x, y ^ x),
    }
}

/// Exact group law, every special case (`P + O`, `P + (−P)`, doubling,
/// the 2-torsion point `(0, 1)`), one inversion per call.
pub fn add(gf: &Gf2, p: Pt, q: Pt) -> Pt {
    match (p, q) {
        (Pt::Inf, r) | (r, Pt::Inf) => r,
        (Pt::Aff(x1, y1), Pt::Aff(x2, y2)) => {
            if x1 == x2 {
                if y1 ^ y2 == x1 {
                    return Pt::Inf;
                }
                // y1 == y2 and x1 != 0: doubling.
                let lam = x1 ^ gf.mul(y1, gf.inv(x1));
                let x3 = gf.sqr(lam) ^ lam ^ 1;
                let y3 = gf.sqr(x1) ^ gf.mul(lam ^ 1, x3);
                return Pt::Aff(x3, y3);
            }
            let lam = gf.mul(y1 ^ y2, gf.inv(x1 ^ x2));
            let x3 = gf.sqr(lam) ^ lam ^ x1 ^ x2 ^ 1;
            let y3 = gf.mul(lam, x1 ^ x3) ^ x3 ^ y1;
            Pt::Aff(x3, y3)
        }
    }
}

/// Whether `(x, y)` is on `y² + xy = x³ + x² + 1`.
pub fn on_curve(gf: &Gf2, x: u64, y: u64) -> bool {
    let x2 = gf.sqr(x);
    gf.sqr(y) ^ gf.mul(x, y) == gf.mul(x2, x) ^ x2 ^ 1
}

/// The rational `y` with `(x, y)` on the curve, the smaller of the two, or
/// `None` when `x` is the abscissa of no rational point.
pub fn lift(gf: &Gf2, solver: &ArtinSchreier, x: u64) -> Option<u64> {
    if x == 0 {
        return Some(1);
    }
    // y = x·z with z² + z = x + 1 + 1/x².
    let c = x ^ 1 ^ gf.sqr(gf.inv(x));
    let z = solver.solve(c)?;
    let y = gf.mul(x, z);
    let y = y.min(y ^ x);
    debug_assert!(on_curve(gf, x, y));
    Some(y)
}

/// The abscissae of `V ∩ {rational}` with a chosen `y` each, in increasing
/// order of `x`.
pub fn base_points(gf: &Gf2, solver: &ArtinSchreier, l: u32) -> Vec<(u64, u64)> {
    (0..(1u64 << l))
        .filter_map(|x| lift(gf, solver, x).map(|y| (x, y)))
        .collect()
}

/// A target for which the group-law formulation can disagree with the
/// baseline: `x_R = 0` (R is 2-torsion), `x_R ∈ V`, or `1/x_R ∈ V`.
pub fn degenerate_target(gf: &Gf2, l: u32, xr: u64) -> bool {
    xr == 0 || xr < (1u64 << l) || gf.inv(xr) < (1u64 << l)
}

/// The four candidate third abscissae for the pair `(P₁, P₂)`:
/// `x(P₁+P₂+R), x(P₁+P₂−R), x(P₁−P₂+R), x(P₁−P₂−R)`, `None` where the sum
/// is the point at infinity.  Exact reference, slow.
pub fn thirds(gf: &Gf2, p1: (u64, u64), p2: (u64, u64), r: (u64, u64)) -> [Option<u64>; 4] {
    let (p1, p2) = (Pt::Aff(p1.0, p1.1), Pt::Aff(p2.0, p2.1));
    let r = Pt::Aff(r.0, r.1);
    let mut out = [None; 4];
    for (k, z) in [add(gf, p1, p2), add(gf, p1, neg(p2))]
        .into_iter()
        .enumerate()
    {
        for (m, rr) in [r, neg(r)].into_iter().enumerate() {
            if let Pt::Aff(x, _) = add(gf, z, rr) {
                out[2 * k + m] = Some(x);
            }
        }
    }
    out
}

/// Every sorted triple `(a ≤ b ≤ c)` of base abscissae with
/// `±P_a ± P_b ± P_c = ±R`, found pair by pair with the exact reference.
/// For a non-[`degenerate_target`] this is the baseline's root set.
pub fn roots_reference(gf: &Gf2, base: &[(u64, u64)], l: u32, r: (u64, u64)) -> Vec<[u64; 3]> {
    let span = 1u64 << l;
    let mut out = Vec::new();
    for i in 0..base.len() {
        for j in i..base.len() {
            for t in thirds(gf, base[i], base[j], r).into_iter().flatten() {
                if t < span {
                    let mut tr = [base[i].0, base[j].0, t];
                    tr.sort_unstable();
                    out.push(tr);
                }
            }
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}
