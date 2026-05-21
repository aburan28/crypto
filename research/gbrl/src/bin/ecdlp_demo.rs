// Toy ECDLP solver via Semaev-style index calculus.
//
// Pipeline:
//   1. Compute |E|, factor it, pick prime subgroup order ℓ.
//   2. Find a generator P of order ℓ.
//   3. Build a factor base F = {F_0, ..., F_{k-1}} of curve points.
//   4. Sample random a in [1, ℓ), compute R = aP, try to decompose R = s_i F_i + s_j F_j
//      via Semaev S_3 + Gröbner basis + group-law filter.
//   5. Each successful decomp gives a linear relation
//        a ≡ s_i · log(F_i) + s_j · log(F_j) (mod ℓ).
//   6. Collect ≥ |F| relations, solve the linear system mod ℓ to recover all log(F_i).
//   7. Verify by brute-force discrete log on each F_i.
//
//     cargo run --release --bin ecdlp_demo

use gbrl::buchberger::BuchbergerState;
use gbrl::curve::{add, filter_decomp_signed, group_order, point_at_x, scalar_mul, Point};
use gbrl::field::{Fp, P};
use gbrl::monomial::Monomial;
use gbrl::poly::{Poly, Term};
use gbrl::semaev::{decomposition_system_2, Curve};
use gbrl::strategy::{Strategy, SugarStrategy};

fn factor_base_poly(values: &[Fp]) -> Poly {
    let mut acc = Poly::from_terms(
        vec![Term { coef: Fp::one(), mono: Monomial::one(1) }],
        1,
    );
    for &v in values {
        let factor = Poly::from_terms(
            vec![
                Term { coef:  Fp::one(), mono: Monomial { exp: vec![1] } },
                Term { coef: -v,         mono: Monomial { exp: vec![0] } },
            ],
            1,
        );
        acc = acc.mul(&factor);
    }
    acc
}

fn trial_factor(n: u64) -> Vec<(u64, u32)> {
    let mut out = vec![];
    let mut n = n;
    let mut d = 2u64;
    while d * d <= n {
        let mut e = 0;
        while n % d == 0 {
            n /= d;
            e += 1;
        }
        if e > 0 { out.push((d, e)); }
        d += 1;
    }
    if n > 1 { out.push((n, 1)); }
    out
}

// Attempt to decompose R into the factor base. Returns canonical-factor-base
// indices (i, j) and the signs (s_i, s_j) such that s_i·F_i + s_j·F_j = R.
fn try_decompose(
    curve: &Curve, x_r: Fp, fb_points: &[Point], fb_poly: &Poly, v_values: &[Fp], r: Point,
) -> Option<(usize, usize, i64, i64)> {
    let sys = decomposition_system_2(curve, x_r, fb_poly);
    let mut state = BuchbergerState::new(sys);
    let mut strategy: Box<dyn Strategy> = Box::new(SugarStrategy);
    while !state.is_done() {
        let idx = strategy.select(&state);
        state.step(idx);
    }
    for (i, x1) in v_values.iter().enumerate() {
        for (j, x2) in v_values.iter().enumerate() {
            let zero = state.basis.iter().all(|p| p.eval(&[*x1, *x2]).is_zero());
            if !zero { continue; }
            if let Some((s1, s2)) = filter_decomp_signed(curve, fb_points[i], fb_points[j], r) {
                return Some((i, j, s1, s2));
            }
        }
    }
    None
}

// Brute-force discrete log of `q` wrt base `p` of order `ell`. O(ell).
fn brute_log(curve: &Curve, p: Point, ell: u64, q: Point) -> Option<u64> {
    let mut acc = Point::Infinity;
    for k in 0..ell {
        if acc == q { return Some(k); }
        acc = add(curve, acc, p);
    }
    None
}

// Modular inverse via extended Euclidean. Returns None if not coprime.
fn mod_inv(a: i64, m: i64) -> Option<i64> {
    let a = a.rem_euclid(m);
    if a == 0 { return None; }
    let (mut old_r, mut r) = (a, m);
    let (mut old_s, mut s) = (1i64, 0i64);
    while r != 0 {
        let q = old_r / r;
        let (nr, ns) = (old_r - q * r, old_s - q * s);
        old_r = r; r = nr;
        old_s = s; s = ns;
    }
    if old_r != 1 { return None; }
    Some(old_s.rem_euclid(m))
}

// Gaussian elimination over Z/ell Z on M x = b. Mutates inputs in-place.
// Returns Some(x) if a unique solution exists; None if singular / underdetermined.
fn solve_mod(matrix: &mut Vec<Vec<i64>>, b: &mut Vec<i64>, ell: i64) -> Option<Vec<i64>> {
    let n_rows = matrix.len();
    let n_cols = matrix[0].len();
    assert_eq!(b.len(), n_rows);

    let mut row = 0;
    for col in 0..n_cols {
        if row >= n_rows { break; }
        let mut pivot = None;
        for r in row..n_rows {
            if matrix[r][col].rem_euclid(ell) != 0 {
                pivot = Some(r); break;
            }
        }
        let pr = match pivot { Some(r) => r, None => continue };
        matrix.swap(row, pr);
        b.swap(row, pr);

        let pv = matrix[row][col].rem_euclid(ell);
        let pv_inv = mod_inv(pv, ell)?;
        for c in 0..n_cols {
            matrix[row][c] = (matrix[row][c] * pv_inv).rem_euclid(ell);
        }
        b[row] = (b[row] * pv_inv).rem_euclid(ell);

        for r in 0..n_rows {
            if r == row { continue; }
            let factor = matrix[r][col].rem_euclid(ell);
            if factor == 0 { continue; }
            for c in 0..n_cols {
                matrix[r][c] = (matrix[r][c] - factor * matrix[row][c]).rem_euclid(ell);
            }
            b[r] = (b[r] - factor * b[row]).rem_euclid(ell);
        }
        row += 1;
    }
    if row < n_cols { return None; }
    let mut x = vec![0i64; n_cols];
    for c in 0..n_cols { x[c] = b[c]; }
    Some(x)
}

fn main() {
    let curve = Curve { a: Fp::new(1), b: Fp::new(1) };
    println!("curve  y² = x³ + x + 1   over F_{}", P);

    let order = group_order(&curve);
    println!("|E|    = {}", order);
    let factors = trial_factor(order);
    println!("|E| = {}", factors.iter()
        .map(|(p, e)| if *e == 1 { format!("{}", p) } else { format!("{}^{}", p, e) })
        .collect::<Vec<_>>().join(" · "));

    // Pick the largest prime factor with order between 30 and 5000 (small enough
    // for brute_log to verify, big enough to be a real subgroup).
    let ell = factors.iter()
        .filter(|(p, e)| *e == 1 && *p > 30 && *p < 5000)
        .map(|(p, _)| *p)
        .max()
        .or_else(|| factors.iter().map(|(p, _)| *p).max())
        .expect("at least one factor");
    let cofactor = order / ell;
    println!("chose ℓ = {} (cofactor {})", ell, cofactor);

    // Find generator P of order ell.
    let mut p_gen = Point::Infinity;
    for x_try in 0..2000 {
        let x = Fp::new(x_try as i64);
        if let Some(pt) = point_at_x(&curve, x) {
            let g = scalar_mul(&curve, pt, cofactor);
            if g != Point::Infinity && scalar_mul(&curve, g, ell) == Point::Infinity {
                p_gen = g;
                break;
            }
        }
    }
    assert!(p_gen != Point::Infinity, "no generator found");
    println!("generator P = {:?}", p_gen);

    // Build factor base IN the subgroup <P>. For each candidate raw curve point,
    // project onto <P> via cofactor multiplication: F'_i = cofactor · raw_i. Then
    // F'_i ∈ <P> (since its order divides ℓ); the projection is the identity iff
    // raw_i was already in <P>, and Infinity iff raw_i's order is coprime to ℓ.
    // Use the surviving projected points as the factor base.
    const FB_SIZE: usize = 12;
    let mut fb_points: Vec<Point> = Vec::new();
    for x_try in 0..30_000u64 {
        if fb_points.len() == FB_SIZE { break; }
        let x = Fp::new(x_try as i64);
        if let Some(raw) = point_at_x(&curve, x) {
            let proj = scalar_mul(&curve, raw, cofactor);
            if proj == Point::Infinity { continue; }
            // Dedup by x-coord (different raws can project to the same point).
            let proj_x = match proj { Point::Affine(x, _) => x, _ => continue };
            if fb_points.iter().any(|p| matches!(p, Point::Affine(x, _) if *x == proj_x)) {
                continue;
            }
            fb_points.push(proj);
        }
    }
    assert_eq!(fb_points.len(), FB_SIZE, "couldn't build factor base in subgroup");
    let v_values: Vec<Fp> = fb_points.iter()
        .map(|p| match p { Point::Affine(x, _) => *x, _ => unreachable!() })
        .collect();
    println!("factor base ({} pts) at x = {:?}", fb_points.len(),
        v_values.iter().map(|f| f.0).collect::<Vec<_>>());
    let fb_poly = factor_base_poly(&v_values);

    // Pre-compute true logs to validate later (cheating, but only for the check).
    println!("\ncomputing true logs of factor base (brute force, ℓ = {})...", ell);
    let true_logs: Vec<Option<u64>> = fb_points.iter()
        .map(|fp| brute_log(&curve, p_gen, ell, *fp))
        .collect();
    let n_in_subgroup = true_logs.iter().filter(|t| t.is_some()).count();
    println!("  {} of {} factor base points lie in <P> (subgroup of order {})",
        n_in_subgroup, FB_SIZE, ell);
    for (i, tl) in true_logs.iter().enumerate() {
        println!("  F_{:>2} = {:?}  true log = {:?}", i, fb_points[i], tl);
    }

    // Collect relations from random a in [1, ell). Need substantially more than
    // FB_SIZE so the relation matrix has full rank — some basis vars take many
    // tries to appear in a successful decomposition.
    let need = FB_SIZE * 4;
    let mut relations: Vec<(u64, Vec<i64>)> = vec![];
    let mut tried = 0usize;
    let mut a: u64 = 1;
    while relations.len() < need && tried < 20_000 {
        tried += 1;
        let r = scalar_mul(&curve, p_gen, a);
        if r == Point::Infinity { a += 1; continue; }
        let x_r = match r { Point::Affine(x, _) => x, _ => unreachable!() };
        if let Some((i, j, si, sj)) = try_decompose(&curve, x_r, &fb_points, &fb_poly, &v_values, r) {
            let mut row = vec![0i64; FB_SIZE];
            row[i] += si;
            row[j] += sj;
            // If signs cancel for i == j, this is a degenerate row; skip.
            if row.iter().all(|&v| v == 0) { a += 1; continue; }
            relations.push((a, row));
            if relations.len() <= 16 || relations.len() % 5 == 0 {
                println!("  relation #{:>2}: a = {:>4}  → {:+} F_{} {:+} F_{}",
                    relations.len(), a, si, i, sj, j);
            }
        }
        a += 1;
    }
    println!("collected {} relations in {} trials", relations.len(), tried);
    if relations.len() < FB_SIZE {
        println!("not enough relations — increase trial budget or factor base.");
        return;
    }

    // Build and solve M x = b mod ell.
    let mut m: Vec<Vec<i64>> = relations.iter().map(|(_, row)| row.clone()).collect();
    let mut b: Vec<i64> = relations.iter().map(|(a, _)| *a as i64).collect();
    println!("\nsolving {}x{} linear system mod {}...", m.len(), FB_SIZE, ell);
    let xs = match solve_mod(&mut m, &mut b, ell as i64) {
        Some(x) => x,
        None => {
            println!("linear system singular — collect more relations.");
            return;
        }
    };

    println!("\nrecovered logs (mod {}):", ell);
    let mut correct = 0;
    let mut subgroup_count = 0;
    for (i, x) in xs.iter().enumerate() {
        let tl = true_logs[i];
        let in_subgroup = tl.is_some();
        if in_subgroup { subgroup_count += 1; }
        let ok = match tl {
            Some(t) => (t as i64) == *x,
            None => false, // F_i not in <P>, recovered "log" is meaningless
        };
        if ok { correct += 1; }
        let marker = match tl {
            Some(_) => if ok { "✓" } else { "✗" },
            None => "—",
        };
        println!("  log(F_{:>2}) = {:>5}   true = {:?}   {}", i, x, tl, marker);
    }
    println!("\n{}/{} factor base points lie in <P>, of which {}/{} logs recovered correctly.",
        subgroup_count, FB_SIZE, correct, subgroup_count);
}
