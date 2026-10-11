//! Exact decomposition census for Trimoska's thesis Table 7.7.
//!
//! Curve `E_1 : y² + xy = x³ + x² + 1` over `F_{2^24}` (the thesis' Koblitz
//! curve), factor base `F_ℓ = {P : x(P) ∈ V_ℓ}`, `V_ℓ = ⟨1, z, …, z^{ℓ−1}⟩`
//! in the polynomial basis, `m = 3`.  For every unordered triple of
//! factor-base points the sum `R` is formed and `x(R)` counted, so for every
//! abscissa `x` of the field we know exactly how many x-triples the x-only
//! `S₄` system admits.  The thesis' `P` is the fraction of random targets
//! whose x-only system is satisfiable; it is compared with
//! `P_approx = 2^{mℓ}/(m!·2^n)` and with the Poisson prediction
//! `1 − exp(−λ)`, `λ = C(|F|, 3)/#E`, on the curve and on its twist.
//!
//! No search is involved: this is a census, `Θ(|F|³)` group additions.
//! The field is the repository's `u64` polynomial basis (the first
//! irreducible of degree 24 in integer order); the thesis generator's
//! modulus is not vendored, and the census is a statement about counts
//! that do not depend on the modulus beyond per-curve fluctuation.

use std::collections::HashMap;

fn clmul_mod(mut a: u64, mut b: u64, irr: u64, n: u32) -> u64 {
    let mut r = 0u64;
    let top = 1u64 << n;
    while b != 0 {
        if b & 1 == 1 {
            r ^= a;
        }
        b >>= 1;
        a <<= 1;
        if a & top != 0 {
            a ^= irr;
        }
    }
    r
}

fn is_irreducible(f: u64, n: u32) -> bool {
    let powx = |k: u32| -> u64 {
        let mut h = 2u64;
        for _ in 0..k {
            h = clmul_mod(h, h, f, n);
        }
        h
    };
    if powx(n) != 2 {
        return false;
    }
    let mut primes = Vec::new();
    let mut m = n;
    let mut p = 2;
    while p * p <= m {
        if m % p == 0 {
            primes.push(p);
            while m % p == 0 {
                m /= p;
            }
        }
        p += 1;
    }
    if m > 1 {
        primes.push(m);
    }
    for p in primes {
        let h = powx(n / p) ^ 2;
        // gcd(h, f) must be 1
        let (mut a, mut b) = (f, h);
        while b != 0 {
            // a mod b
            let db = 63 - b.leading_zeros();
            while a != 0 && (63 - a.leading_zeros()) >= db {
                a ^= b << ((63 - a.leading_zeros()) - db);
            }
            std::mem::swap(&mut a, &mut b);
        }
        if a != 1 {
            return false;
        }
    }
    true
}

fn find_irreducible(n: u32) -> u64 {
    let top = 1u64 << n;
    for low in 1u64..top {
        let f = top | low;
        if f & 1 == 1 && is_irreducible(f, n) {
            return f;
        }
    }
    unreachable!()
}

struct Field {
    n: u32,
    irr: u64,
    /// Row-reduced solver for z² + z = c: `solve_tab[i]` applied to c gives z.
    lin: LinSolver,
}

/// Solves the F2-linear system L(z) = c, L(z) = z² + z, by precomputed elimination.
struct LinSolver {
    n: usize,
    /// For each pivot column `col`, the mask (over c-bits) whose parity gives z_col.
    masks: Vec<(usize, u64)>,
    /// Consistency mask: c is in the image iff parity(c & kernel_check) == 0.
    kernel_check: u64,
}

impl LinSolver {
    fn new(n: usize, l: impl Fn(u64) -> u64) -> Self {
        // Build rows: row j = (bit j of L(e_i) for all i) | identity part.
        let mut rows: Vec<(u64, u64)> = (0..n)
            .map(|j| {
                let mut row = 0u64;
                for i in 0..n {
                    if (l(1u64 << i) >> j) & 1 == 1 {
                        row |= 1u64 << i;
                    }
                }
                (row, 1u64 << j)
            })
            .collect();
        let mut pivots: Vec<(usize, usize)> = Vec::new();
        let mut r = 0;
        for col in 0..n {
            if let Some(p) = (r..n).find(|&i| (rows[i].0 >> col) & 1 == 1) {
                rows.swap(r, p);
                let piv = rows[r];
                for i in 0..n {
                    if i != r && (rows[i].0 >> col) & 1 == 1 {
                        rows[i].0 ^= piv.0;
                        rows[i].1 ^= piv.1;
                    }
                }
                pivots.push((col, r));
                r += 1;
            }
        }
        // Read the transforms only after every column has been cleared.
        let masks: Vec<(usize, u64)> = pivots
            .iter()
            .map(|&(col, row)| (col, rows[row].1))
            .collect();
        assert_eq!(r, n - 1, "z^2+z has a kernel of dimension exactly 1");
        let kernel_check = rows[n - 1].1;
        LinSolver {
            n,
            masks,
            kernel_check,
        }
    }
    fn solve(&self, c: u64) -> Option<u64> {
        if (c & self.kernel_check).count_ones() & 1 == 1 {
            return None;
        }
        let mut z = 0u64;
        for &(col, mask) in &self.masks {
            if (c & mask).count_ones() & 1 == 1 {
                z |= 1u64 << col;
            }
        }
        debug_assert!(self.n <= 64);
        Some(z)
    }
}

impl Field {
    fn new(n: u32) -> Field {
        let irr = find_irreducible(n);
        let lin = LinSolver::new(n as usize, |z| clmul_mod(z, z, irr, n) ^ z);
        Field { n, irr, lin }
    }
    #[inline]
    fn mul(&self, a: u64, b: u64) -> u64 {
        clmul_mod(a, b, self.irr, self.n)
    }
    #[inline]
    fn sqr(&self, a: u64) -> u64 {
        self.mul(a, a)
    }
    fn inv(&self, a: u64) -> u64 {
        assert!(a != 0);
        let mut r = 1u64;
        let mut base = a;
        let mut e = (1u64 << self.n) - 2;
        while e > 0 {
            if e & 1 == 1 {
                r = self.mul(r, base);
            }
            base = self.sqr(base);
            e >>= 1;
        }
        r
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Pt {
    Inf,
    Aff(u64, u64),
}

/// `y² + xy = x³ + a x² + b`.
struct Curve<'a> {
    f: &'a Field,
    a: u64,
    b: u64,
}

impl<'a> Curve<'a> {
    fn lift(&self, x: u64) -> Option<[Pt; 2]> {
        let f = self.f;
        if x == 0 {
            // y² = b
            let mut y = self.b;
            for _ in 0..f.n - 1 {
                y = f.sqr(y);
            }
            return Some([Pt::Aff(0, y), Pt::Aff(0, y)]);
        }
        // (y/x)² + (y/x) = x + a + b/x²
        let inv = f.inv(x);
        let c = x ^ self.a ^ f.mul(self.b, f.sqr(inv));
        let z = f.lin.solve(c)?;
        let y = f.mul(x, z);
        Some([Pt::Aff(x, y), Pt::Aff(x, y ^ x)])
    }
    fn on_curve(&self, p: Pt) -> bool {
        match p {
            Pt::Inf => true,
            Pt::Aff(x, y) => {
                let f = self.f;
                f.sqr(y) ^ f.mul(x, y) == f.mul(f.sqr(x), x) ^ f.mul(self.a, f.sqr(x)) ^ self.b
            }
        }
    }
    fn add(&self, p: Pt, q: Pt) -> Pt {
        let f = self.f;
        match (p, q) {
            (Pt::Inf, _) => q,
            (_, Pt::Inf) => p,
            (Pt::Aff(x1, y1), Pt::Aff(x2, y2)) => {
                if x1 == x2 {
                    if y1 != y2 || x1 == 0 {
                        return Pt::Inf;
                    }
                    let lam = x1 ^ f.mul(y1, f.inv(x1));
                    let x3 = f.sqr(lam) ^ lam ^ self.a;
                    let y3 = f.sqr(x1) ^ f.mul(lam ^ 1, x3);
                    Pt::Aff(x3, y3)
                } else {
                    let lam = f.mul(y1 ^ y2, f.inv(x1 ^ x2));
                    let x3 = f.sqr(lam) ^ lam ^ x1 ^ x2 ^ self.a;
                    let y3 = f.mul(lam, x1 ^ x3) ^ x3 ^ y1;
                    Pt::Aff(x3, y3)
                }
            }
        }
    }
}

fn choose3(k: f64) -> f64 {
    k * (k - 1.0) * (k - 2.0) / 6.0
}

fn census(f: &Field, curve: &Curve, l: u32, label: &str, n_abscissae: usize, order: u64) -> String {
    let n = f.n;
    // Factor base: all points with x in V_l = span(1, z, ..., z^{l-1}) = {0 .. 2^l - 1}.
    let mut base: Vec<Pt> = Vec::new();
    let mut base_x = 0usize;
    for x in 0..(1u64 << l) {
        if let Some([p, q]) = curve.lift(x) {
            base_x += 1;
            base.push(p);
            if q != p {
                base.push(q);
            }
        }
    }
    let size = base.len();
    // Pair sums first, then add the third point: Θ(|F|³/6) additions.
    let mut counts: Vec<u8> = vec![0; 1usize << n];
    let mut counts_distinct_x: Vec<u8> = vec![0; 1usize << n];
    let mut triples = 0u64;
    for i in 0..size {
        for j in (i + 1)..size {
            let pij = curve.add(base[i], base[j]);
            for k in (j + 1)..size {
                let r = curve.add(pij, base[k]);
                assert!(curve.on_curve(r), "addition left the curve");
                triples += 1;
                if let Pt::Aff(x, _) = r {
                    let c = &mut counts[x as usize];
                    *c = c.saturating_add(1);
                    let (Pt::Aff(xi, _), Pt::Aff(xj, _), Pt::Aff(xk, _)) =
                        (base[i], base[j], base[k])
                    else {
                        unreachable!()
                    };
                    if xi != xj && xj != xk && xi != xk {
                        let d = &mut counts_distinct_x[x as usize];
                        *d = d.saturating_add(1);
                    }
                }
            }
        }
    }
    let total: u64 = counts.iter().map(|&c| c as u64).sum();
    // Classify abscissae on this curve: those that lift.
    let mut lifting = 0u64;
    let mut hit = 0u64;
    let mut hit_distinct = 0u64;
    let mut sum = 0u64;
    let mut hist: HashMap<u8, u64> = HashMap::new();
    for x in 0..(1u64 << n) {
        if curve.lift(x).is_some() {
            lifting += 1;
            let c = counts[x as usize];
            sum += c as u64;
            *hist.entry(c.min(9)).or_insert(0) += 1;
            if c > 0 {
                hit += 1;
            }
            if counts_distinct_x[x as usize] > 0 {
                hit_distinct += 1;
            }
        }
    }
    let p_exact = hit as f64 / lifting as f64;
    let p_exact_distinct = hit_distinct as f64 / lifting as f64;
    let lambda_measured = sum as f64 / lifting as f64;
    let lambda_model = choose3(size as f64) / order as f64;
    let p_poisson = 1.0 - (-lambda_model).exp();
    let p_approx = 2f64.powi((3 * l) as i32) / (6.0 * 2f64.powi(n as i32));
    let hist_s: String = (0..=9u8)
        .map(|k| format!("{}:{}", k, hist.get(&k).copied().unwrap_or(0)))
        .collect::<Vec<_>>()
        .join(" ");
    format!(
        "{label} l={l} sum_counts_all_x={total} |V_x lifting|={base_x} |F|={size} triples={triples} order={order} lifting_x={lifting}/{n_abscissae} \
         P_exact={p_exact:.5} P_exact_distinct_x={p_exact_distinct:.5} lambda_measured={lambda_measured:.5} \
         lambda_model=C(|F|,3)/#E={lambda_model:.5} P_poisson={p_poisson:.5} P_approx=2^(3l)/(3!2^n)={p_approx:.5} \
         ratio_exact/approx={:.3} ratio_exact/poisson={:.3} hist[{hist_s}]",
        p_exact / p_approx,
        p_exact / p_poisson
    )
}

fn curve_order(f: &Field, curve: &Curve) -> u64 {
    // 1 (infinity) + points over all x.
    let mut count = 1u64;
    for x in 0..(1u64 << f.n) {
        if let Some([p, q]) = curve.lift(x) {
            assert!(
                curve.on_curve(p) && curve.on_curve(q),
                "lift produced an off-curve point at x={x:#x}"
            );
            count += if p == q { 1 } else { 2 };
        }
    }
    count
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n: u32 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(24);
    let ls: Vec<u32> = args
        .get(2)
        .map(|s| s.split(',').map(|t| t.parse().unwrap()).collect())
        .unwrap_or_else(|| vec![4, 5, 6, 7, 8]);
    let f = Field::new(n);
    println!("field n={n} irreducible=0x{:x}", f.irr);
    // The thesis curve: a = 1, b = 1.  Its quadratic twist: a' = a + delta, Tr(delta) = 1.
    let curve = Curve { f: &f, a: 1, b: 1 };
    let delta = (1..(1u64 << n))
        .find(|&d| {
            let mut s = 0;
            let mut v = d;
            for _ in 0..n {
                s ^= v;
                v = f.sqr(v);
            }
            s == 1
        })
        .unwrap();
    let twist = Curve {
        f: &f,
        a: 1 ^ delta,
        b: 1,
    };
    let order = curve_order(&f, &curve);
    let twist_order = curve_order(&f, &twist);
    println!("curve y^2+xy=x^3+x^2+1: #E={order}  twist (a=1+delta, delta=0x{delta:x}): #E'={twist_order}  sum={} (expected 2*(2^n+1)={})", order + twist_order, 2 * ((1u64 << n) + 1));
    let n_abscissae = 1usize << n;
    for &l in &ls {
        println!("{}", census(&f, &curve, l, "curve", n_abscissae, order));
        println!(
            "{}",
            census(&f, &twist, l, "twist", n_abscissae, twist_order)
        );
    }
}
