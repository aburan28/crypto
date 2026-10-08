use isogeny_algos::field::Rng;
use isogeny_algos::int::Int;
use isogeny_algos::quat::*;

fn rand_o0_elt(alg: &Alg, o0: &Lattice, rng: &mut Rng, r: i64) -> Quat {
    let b = o0.basis();
    let mut x = Quat::zero();
    for e in &b {
        let c = rng.below(2 * r as u64 + 1) as i64 - r;
        x = x.add(&e.scale(&Int::from(c), &Int::one()));
    }
    let _ = alg;
    x
}

/// Left O_0-ideal O_0 n + O_0 a.
fn ideal(alg: &Alg, o0: &Lattice, n: &Int, a: &Quat) -> Lattice {
    let mut gens = vec![];
    for b in o0.basis() {
        gens.push(b.scale(n, &Int::one()));
        gens.push(alg.mul(&b, a));
    }
    Lattice::from_gens(&gens)
}

#[test]
fn maximal_order_and_ideals() {
    let mut rng = Rng::new(900);
    for p in [7i64, 11, 103, 1_000_003, 2_147_483_647] {
        let alg = Alg::new(&Int::from(p));
        let o0 = alg.o0();
        assert_eq!(o0.det(), (Int::one(), Int::from(4i64)), "covolume 1/4");
        // ring: products of basis elements stay in O_0
        for a in o0.basis() {
            for b in o0.basis() {
                assert!(o0.contains(&alg.mul(&a, &b)));
            }
        }
        assert!(o0.contains(&Quat::one()));
        // reduced norm is multiplicative and integral on O_0
        for _ in 0..20 {
            let (x, y) = (
                rand_o0_elt(&alg, &o0, &mut rng, 50),
                rand_o0_elt(&alg, &o0, &mut rng, 50),
            );
            let (nx, dx) = alg.nrd(&x);
            let (ny, dy) = alg.nrd(&y);
            assert!(
                dx == Int::one() && dy == Int::one(),
                "integral norms on O_0"
            );
            assert_eq!(alg.nrd(&alg.mul(&x, &y)), (&nx * &ny, Int::one()));
        }
        // ideals of prime norm
        for &n in &[3i64, 5, 7, 13] {
            if n == p {
                continue;
            }
            let ni = Int::from(n);
            let a = loop {
                let a = rand_o0_elt(&alg, &o0, &mut rng, 20);
                let (na, _) = alg.nrd(&a);
                if !a.is_zero()
                    && na.modulo(&ni).is_zero()
                    && !(o0.scale(&ni, &Int::one())).contains(&a)
                {
                    break a;
                }
            };
            let i = ideal(&alg, &o0, &ni, &a);
            assert_eq!(i.norm(&o0), (ni.clone(), Int::one()), "N(I) = {n}, p = {p}");
            assert_eq!(i.left_order(&alg, &o0), o0);
            let ro = i.right_order(&alg, &o0);
            assert_eq!(ro.det(), o0.det(), "right order is maximal");
            for x in ro.basis() {
                for y in ro.basis() {
                    assert!(ro.contains(&alg.mul(&x, &y)));
                }
            }
            // I Ibar = N(I) O_0
            assert_eq!(i.mul(&alg, &i.conj()), o0.scale(&ni, &Int::one()));
            // LLL keeps the lattice
            let rb = i.reduced_basis(&alg);
            let l2 = Lattice::from_gens(
                &rb.iter()
                    .map(|r| Quat::new(r.clone(), i.den.clone()))
                    .collect::<Vec<_>>(),
            );
            assert_eq!(l2, i);
        }
        // units of O_0: norm-1 elements (scaled bound: den^2 * 1)
        let sv = o0.short_vectors(&alg, &(&o0.den * &o0.den), 100);
        let units = sv.iter().filter(|(_, n)| *n == &o0.den * &o0.den).count();
        let expect = if p == 3 { 6 } else { 2 };
        assert_eq!(units, expect, "units up to sign, p = {p}");
    }
}

#[test]
fn brandt_matrix_spectrum_matches_supersingular_graph() {
    use isogeny_algos::quat::brandt::*;
    let mut rng = Rng::new(901);
    for p in [103i64, 431, 1019] {
        let alg = Alg::new(&Int::from(p));
        for ell in [2u64, 3] {
            let cs = class_set(&alg, ell, &mut rng);
            let (js, a) = supersingular_graph(p as u64, ell as usize, &mut rng);
            let h = cs.reps.len();
            assert_eq!(
                h,
                js.len(),
                "class number = number of supersingular j, p = {p}"
            );
            // expected class number: floor(p/12) + {0,1,1,2} for p = 1,5,7,11 mod 12
            let extra = [0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 2][(p % 12) as usize];
            assert_eq!(h as i64, p / 12 + extra);
            for row in &cs.brandt {
                assert_eq!(row.iter().sum::<u32>(), ell as u32 + 1);
            }
            assert_eq!(
                power_traces(&cs.brandt, h),
                power_traces(&a, h),
                "spectra, p = {p}, l = {ell}"
            );
        }
    }
}

/// A left O_0-ideal of norm n (prime or not): O_0 n + O_0 alpha with n | Nrd(alpha), alpha found
/// by a modular square root when n is prime, by sampling otherwise.
fn random_ideal(alg: &Alg, o0: &Lattice, n: &Int, rng: &mut Rng) -> Lattice {
    use isogeny_algos::quat::klpt::ideal_from;
    loop {
        let x = if n.is_probable_prime() {
            let (b, c, d) = (
                Int::from(rng.next() >> 2),
                Int::from(rng.next() >> 2),
                Int::from(rng.next() >> 2),
            );
            let t = -&(&(&b * &b) + &(&alg.p * &(&(&c * &c) + &(&d * &d))));
            let Some(a) = Int::sqrt_mod_prime(&t, n) else {
                continue;
            };
            Quat::new([a, b, c, d], Int::one())
        } else {
            let x = rand_o0_elt(alg, o0, rng, 1000);
            if x.is_zero() || !alg.nrd(&x).0.modulo(n).is_zero() {
                continue;
            }
            x
        };
        let i = ideal_from(alg, o0, n, &x);
        if i.norm(o0) == (n.clone(), Int::one()) {
            return i;
        }
    }
}

#[test]
fn klpt_outputs_equivalent_two_power_norm_ideals() {
    use isogeny_algos::bigint::Big;
    use isogeny_algos::quat::klpt::*;
    let mut rng = Rng::new(902);
    for ps in [
        "2147483647",
        "1152921504606847067",
        "1267650600228229401496703205707",
    ] {
        let p = Int::from_big(&Big::from_dec(ps));
        let alg = Alg::new(&p);
        let o0 = alg.o0();
        let mut big_n = &p + &Int::from(2i64);
        while !big_n.is_probable_prime() {
            big_n = &big_n + &Int::from(2i64);
        }
        for n in [Int::from(3i64.pow(5) * 5), big_n] {
            let i = random_ideal(&alg, &o0, &n, &mut rng);
            let res = klpt(&alg, &i, 2, &mut rng).expect("KLPT");
            assert_eq!(
                res.j.norm(&o0),
                (Int::from(2i64).pow(res.e), Int::one()),
                "p = {p}"
            );
            assert!(o0.contains_lattice(&res.j), "J in O_0");
            assert_eq!(res.j.left_order(&alg, &o0), o0, "left O_0-ideal");
            assert_eq!(i.rmul(&alg, &res.xi), res.j, "J = I xi");
            assert!(
                !o0.scale(&Int::from(2i64), &Int::one())
                    .contains_lattice(&res.j),
                "primitive"
            );
            let logp = p.to_f64().log2();
            eprintln!(
                "p ~ 2^{logp:.0}, N(I) = {n}: e = {} ({:.2} log2 p)",
                res.e,
                res.e as f64 / logp
            );
        }
    }
}

#[test]
fn deuring_classes_to_supersingular_curves() {
    use isogeny_algos::curve::{jinv, pmul, Pt};
    use isogeny_algos::field::Field;
    use isogeny_algos::quat::brandt::*;
    use isogeny_algos::quat::deuring::Deuring;
    use std::collections::HashMap;
    let mut rng = Rng::new(903);
    // p + 1 = 1260 = 4 * 3^2 * 5 * 7, p - 1 = 2 * 17 * 37
    let p = 1259u64;
    let d = Deuring::new(p, &mut rng);
    assert_eq!(d.t_odd, 315 * 629);
    let alg = &d.alg;
    let cs = class_set(alg, 2, &mut rng);
    let (js, a) = supersingular_graph(p, 2, &mut rng);
    let h = js.len();
    assert_eq!(cs.reps.len(), h);
    // Deuring map on class representatives: a bijection onto the supersingular j-invariants
    let mut map = vec![];
    for (t, i) in cs.reps.iter().enumerate() {
        let j = d
            .ideal_to_j(i, &mut rng)
            .unwrap_or_else(|| panic!("class {t}: no smooth equivalent"));
        map.push(j);
    }
    let jidx: HashMap<(u64, u64), usize> = js.iter().enumerate().map(|(k, &j)| (j, k)).collect();
    let perm: Vec<usize> = map
        .iter()
        .map(|j| *jidx.get(j).expect("Deuring image is supersingular"))
        .collect();
    let mut sorted = perm.clone();
    sorted.sort();
    sorted.dedup();
    assert_eq!(sorted.len(), h, "bijection");
    assert_eq!(map[0], (1728 % p, 0), "O_0 -> E_0");
    // Brandt matrix = Phi_2 multiplicity matrix under the bijection
    for r in 0..h {
        for c in 0..h {
            assert_eq!(cs.brandt[r][c], a[perm[r]][perm[c]], "B(2)[{r}][{c}]");
        }
    }
    // kernel -> ideal -> kernel round trip on E_0 (over F_{p^4}, so p - 1 torsion too)
    let f2 = d.f2.clone();
    for n in [3u64, 5, 7, 9, 15, 35, 63, 315, 17, 37, 629, 3 * 17 * 37] {
        for _ in 0..3 {
            // random point of order n
            let k = loop {
                let r = isogeny_algos::curve::random_point_f(&f2, &d.e0, &mut rng);
                let k = pmul(&f2, &d.e0, &r, (p as u128 * p as u128 - 1) / n as u128);
                let ok = isogeny_algos::field::factor_u64(n)
                    .iter()
                    .all(|&(l, _)| pmul(&f2, &d.e0, &k, (n / l) as u128) != Pt::Inf);
                if ok {
                    break k;
                }
            };
            let id = d.ideal_of_kernel(&k, n, &mut rng).expect("ideal of kernel");
            assert_eq!(id.norm(&d.o0), (Int::from(n), Int::one()));
            let ks = d.kernel_of_ideal(&id).expect("kernel of ideal");
            let cod1 = d.isogeny_from_kernels(&ks);
            // isogeny straight from K
            let mut parts = vec![];
            for (l, e) in isogeny_algos::field::factor_u64(n) {
                parts.push((l, e, pmul(&f2, &d.e0, &k, (n / l.pow(e)) as u128)));
            }
            let cod2 = d.isogeny_from_kernels(&parts);
            assert_eq!(jinv(&f2, &cod1), jinv(&f2, &cod2), "n = {n}");
            let _ = f2.one();
        }
    }
}
