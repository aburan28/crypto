//! Ideal classes of O_0, Brandt matrices, and the supersingular l-isogeny graph on the j-line over
//! F_{p^2} (Mestre's "méthode des graphes"). The Brandt matrix B(l)[i][j] counts the l + 1
//! sub-ideals of norm l N(I_i) of the class representative I_i that fall in class j; the j-graph
//! matrix A(l)[i][j] is the multiplicity of j_j as a root of Phi_l(j_i, Y). Deuring's
//! correspondence makes them equal under the class -> j-invariant bijection.
use super::*;
use crate::field::{Field, Rng, Zp, Zp2};
use crate::find::modpoly::Phi;
use crate::poly;
use std::collections::HashMap;

/// beta with J = I beta if the left ideals I, J (of the same left order) are equivalent: an element
/// x of Ibar J of norm N(I) N(J) gives beta = x / N(I).
pub fn equivalent(alg: &Alg, o: &Lattice, i: &Lattice, j: &Lattice) -> Option<Quat> {
    let (ni, di) = i.norm(o);
    let (nj, dj) = j.norm(o);
    let prod = i.conj().mul(alg, j);
    // target Nrd = N(I) N(J); scaled by den^2 of the lattice
    let target_num = &ni * &nj;
    let target_den = &di * &dj;
    let scaled = (&target_num * &(&prod.den * &prod.den)).div_floor(&target_den);
    for (v, nv) in prod.short_vectors(alg, &scaled, 4) {
        if &nv * &target_den == &target_num * &(&prod.den * &prod.den) {
            let x = Quat::new(v, prod.den.clone());
            return Some(x.scale(&di, &ni));
        }
    }
    None
}

/// Number of units of the order up to sign (|O^x| / 2): elements of reduced norm 1.
pub fn units_mod_sign(alg: &Alg, o: &Lattice) -> usize {
    let d2 = &o.den * &o.den;
    o.short_vectors(alg, &d2, 64)
        .iter()
        .filter(|(_, n)| *n == d2)
        .count()
}

/// A reduced representative of the class of I: J = I deltabar / N(I) for a shortest delta in I.
pub fn reduce_ideal(alg: &Alg, o: &Lattice, i: &Lattice) -> Lattice {
    let rb = i.reduced_basis(alg);
    // shortest of the reduced basis vectors (LLL: within 2^{3/2} of the minimum)
    let best = rb.iter().min_by(|a, b| alg.qf(a).cmp(&alg.qf(b))).unwrap();
    let delta = Quat::new(best.clone(), i.den.clone());
    let (n, d) = i.norm(o);
    i.rmul(alg, &delta.conj()).scale(&d, &n)
}

/// The l + 1 left ideals of norm l of a maximal order o (l prime, l != p): o l + o alpha for
/// alpha in o with l | Nrd(alpha), found by sampling small elements.
pub fn ideals_of_norm(alg: &Alg, o: &Lattice, ell: u64, rng: &mut Rng) -> Vec<Lattice> {
    let l = Int::from(ell);
    let lo = o.scale(&l, &Int::one());
    let basis = o.basis();
    let mut out: Vec<Lattice> = vec![];
    let mut tries = 0;
    while out.len() < ell as usize + 1 {
        tries += 1;
        assert!(tries < 100_000, "ideals_of_norm: sampling failed");
        let mut a = Quat::zero();
        for b in &basis {
            let c = rng.below(2 * ell + 1) as i64 - ell as i64;
            a = a.add(&b.scale(&Int::from(c), &Int::one()));
        }
        let (na, da) = alg.nrd(&a);
        if a.is_zero() || da != Int::one() || !na.modulo(&l).is_zero() || lo.contains(&a) {
            continue;
        }
        let mut gens = vec![];
        for b in &basis {
            gens.push(b.scale(&l, &Int::one()));
            gens.push(alg.mul(b, &a));
        }
        let id = Lattice::from_gens(&gens);
        if id.norm(o) == (l.clone(), Int::one()) && !out.contains(&id) {
            out.push(id);
        }
    }
    out
}

pub struct ClassSet {
    /// left O_0-ideal representatives (reduced)
    pub reps: Vec<Lattice>,
    /// |O_R(I)^x| / 2 for each class
    pub units: Vec<usize>,
    /// Brandt matrix B(l): row i = classes of the l + 1 sub-ideals of reps[i]
    pub brandt: Vec<Vec<u32>>,
}

/// Theta-series prefix of a maximal order: the number of elements of each reduced norm up to
/// 4 sqrt(p) (an isomorphism invariant; equivalent left O_0-ideals have conjugate, hence
/// isomorphic, right orders). Used to bucket classes so the exact equivalence test runs only
/// inside a bucket.
pub fn theta_key(alg: &Alg, o: &Lattice) -> Vec<u32> {
    let p = alg.p.to_f64();
    let nmax = (4.0 * p.sqrt()).ceil() as i64 + 4;
    let d2 = &o.den * &o.den;
    let mut counts = vec![0u32; nmax as usize + 1];
    for (_, nv) in o.short_vectors(alg, &(&d2 * &Int::from(nmax)), 1_000_000) {
        let (q, r) = nv.divrem_trunc(&d2);
        if r.is_zero() {
            if let Some(k) = q.to_i128() {
                if (k as usize) < counts.len() {
                    counts[k as usize] += 1;
                }
            }
        }
    }
    counts
}

/// Classes of left O_0-ideals by breadth-first search in the l-neighbour graph, with the Brandt
/// matrix B(l). Classes are bucketed by `theta_key` of the right order, and new ideals are tested
/// for equivalence only against representatives in their bucket. Checks Eichler's mass formula
/// sum 1/|O_R(I)^x| = (p-1)/24 at the end.
pub fn class_set(alg: &Alg, ell: u64, rng: &mut Rng) -> ClassSet {
    let o0 = alg.o0();
    let mut reps = vec![o0.clone()];
    let mut buckets: HashMap<Vec<u32>, Vec<usize>> = HashMap::new();
    buckets.entry(theta_key(alg, &o0)).or_default().push(0);
    let mut rows: Vec<Vec<usize>> = vec![];
    let mut k = 0;
    while k < reps.len() {
        let i = reps[k].clone();
        let ro = i.right_order(alg, &o0);
        let mut row = vec![];
        for l_id in ideals_of_norm(alg, &ro, ell, rng) {
            let j = i.mul(alg, &l_id);
            let jr = reduce_ideal(alg, &o0, &j);
            let key = theta_key(alg, &jr.right_order(alg, &o0));
            let bucket = buckets.entry(key).or_default();
            let idx = bucket
                .iter()
                .copied()
                .find(|&t| equivalent(alg, &o0, &reps[t], &jr).is_some());
            let idx = match idx {
                Some(t) => t,
                None => {
                    reps.push(jr);
                    bucket.push(reps.len() - 1);
                    reps.len() - 1
                }
            };
            row.push(idx);
        }
        rows.push(row);
        k += 1;
    }
    let h = reps.len();
    let brandt = rows
        .iter()
        .map(|r| {
            let mut v = vec![0u32; h];
            for &t in r {
                v[t] += 1;
            }
            v
        })
        .collect();
    let units: Vec<usize> = reps
        .iter()
        .map(|i| units_mod_sign(alg, &i.right_order(alg, &o0)))
        .collect();
    // Eichler: sum 1/|O_i^x| = (p-1)/24 with |O_i^x| = 2 u_i, i.e. 12 sum (W/u_i) = (p-1) W
    // for W = lcm(u_i) | 6
    let w: usize = 6;
    let lhs: usize = units.iter().map(|&u| 12 * (w / u)).sum();
    let p = alg.p.to_i128().expect("p fits i128") as usize;
    assert_eq!(lhs, (p - 1) * w, "Eichler mass formula");
    ClassSet {
        reps,
        units,
        brandt,
    }
}

/// Supersingular j-invariants over F_{p^2} and the Phi_l root-multiplicity matrix, by BFS from
/// j = 1728 (supersingular for p = 3 mod 4).
pub fn supersingular_graph(p: u64, ell: usize, rng: &mut Rng) -> (Vec<(u64, u64)>, Vec<Vec<u32>>) {
    let f2 = Zp2::new(p);
    let phi: Phi<Zp2> = if p > 10_000 {
        Phi::compute(&Zp::new(p), ell).lift(|a| (a, 0))
    } else {
        Phi::via_crt(&f2, ell)
    };
    let start = (1728 % p, 0);
    let mut js = vec![start];
    let mut idx: HashMap<(u64, u64), usize> = HashMap::from([(start, 0)]);
    let mut edges: Vec<Vec<((u64, u64), u32)>> = vec![];
    let mut k = 0;
    while k < js.len() {
        let j = js[k];
        let mut yp = phi.y_poly(&f2, j);
        let mut out = vec![];
        for r in poly::roots(&f2, &yp, rng) {
            let mut m = 0;
            loop {
                let (q, rem) = poly::divrem(&f2, &yp, &vec![f2.neg(r), f2.one()]);
                if !rem.is_empty() {
                    break;
                }
                yp = q;
                m += 1;
            }
            out.push((r, m));
            if let std::collections::hash_map::Entry::Vacant(e) = idx.entry(r) {
                e.insert(js.len());
                js.push(r);
            }
        }
        edges.push(out);
        k += 1;
    }
    let h = js.len();
    let a = edges
        .iter()
        .map(|row| {
            let mut v = vec![0u32; h];
            for &(r, m) in row {
                v[idx[&r]] += m;
            }
            v
        })
        .collect();
    (js, a)
}

/// tr(M^k) mod a 61-bit prime for k = 1..=kmax (equal power sums up to the dimension mean equal
/// characteristic polynomials over a field of characteristic 0 or larger than the dimension).
pub fn power_traces(m: &[Vec<u32>], kmax: usize) -> Vec<u64> {
    const P: u128 = (1u128 << 61) - 1;
    let h = m.len();
    let mut cur: Vec<Vec<u128>> = m
        .iter()
        .map(|r| r.iter().map(|&x| x as u128).collect())
        .collect();
    let mut out = vec![];
    for _ in 0..kmax {
        out.push(((0..h).map(|i| cur[i][i]).sum::<u128>() % P) as u64);
        let mut next = vec![vec![0u128; h]; h];
        for i in 0..h {
            for k in 0..h {
                if cur[i][k] == 0 {
                    continue;
                }
                for j in 0..h {
                    next[i][j] = (next[i][j] + cur[i][k] * m[k][j] as u128) % P;
                }
            }
        }
        cur = next;
    }
    out
}
