//! Isogeny-loop endomorphisms combined with GLS: 4-dimensional GLV-GLS at
//! moderate CM discriminant, measured at the decomposition-lattice level.
//!
//! Follow-up to `RESEARCH_ISOGENY_LOOP_GLV.md` (Koshelev–Sanso, ePrint
//! 2024/1985).  The paper builds a 2-dimensional GLV split over a prime field
//! from a loop of horizontal isogenies.  This example asks whether the same
//! loop endomorphism `alpha` also gives a 4-dimensional split on a GLS curve
//! `E'/F_{p^2}` (the quadratic twist of `E/F_p`), alongside the GLS map `psi`
//! with `psi^2 = -1`.  Everything is integer and lattice arithmetic; no curve
//! is constructed, and that is stated in the output.
//!
//! Three parts:
//!
//! 1. **Cheapest loop per discriminant, certified exactly.**  A loop of
//!    horizontal isogenies is a primitive element `alpha = (A + B sqrt D)/2`
//!    of `O_K`, `B != 0`, with degree `N(alpha) = (A^2 - D B^2)/4`; its prime
//!    ideal factors give one step per prime ideal, each of cost `c(l)`.  In the
//!    paper's projective model `c(l) = 7.5 l`, so `c(l) >= KAPPA log2 l` with
//!    `KAPPA = 7.5*3/log2 3`, and a loop of cost `<= C` has `N <= 2^(C/KAPPA)`.
//!    Enumerating every primitive `alpha` up to that norm certifies the
//!    minimum without a class group (the paper's Section 5 brute force, made
//!    exhaustive by this bound).
//! 2. **GLS instances from integers alone.**  Choose `t, v` with
//!    `p = (t^2 + |D| v^2)/4` prime of `--bits` bits; a curve `E/F_p` with
//!    trace `t` and `End(E) = O_K` exists (CM theory: `p` is a principal prime
//!    of degree one, so it splits completely in the Hilbert class field).  Its
//!    quadratic twist over `F_{p^2}` has `#E'(F_{p^2}) = (p-1)^2 + t^2`; keep
//!    instances where this is `h * r` with `r` prime and `h <= 64`.  On
//!    `E'(F_{p^2})[r]` the GLS map `psi` (the twisted `p`-Frobenius) satisfies
//!    `psi^2 = -1` and `psi^2 - t psi + p = 0`, so its eigenvalue is
//!    `mu = (p-1)/t mod r`.  Transporting `End(E) = O_K` along the twist,
//!    `pi -> mu` fixes the ring map `O_K -> Z/r`: `sqrt D -> s = (2 mu - t)/v`
//!    and `alpha -> lambda = (A + B s)/2`.  The conjugate embedding (`-s`) is
//!    the other eigenvalue; both are reported.  Every identity is asserted:
//!    `mu^2 = -1`, `mu^2 - t mu + p = 0`, `s^2 = D`, `lambda^2 - A lambda + N = 0`.
//! 3. **Decomposition quality.**  LLL (crate `lattice::lll_reduce_hp`, 2048-bit Gram–Schmidt; the f64 variant left 254-bit bases unreduced) on the
//!    lattice `{k : k0 + k1 mu + k2 lambda + k3 lambda mu = 0 mod r}` and on the
//!    GLS-only lattice `{k : k0 + k1 mu = 0}`; Babai rounding of seeded random
//!    scalars; each decomposition's congruence is asserted.  Reported: the max
//!    over samples and the mean of `max_i bits(|k_i|)`, against `log2(r)/4`
//!    and `log2(r)/2`.
//!
//! A derived operation-count model (labelled as such, not measured) prices
//! the two splits with interleaved width-`w` NAF: doublings `c_dbl * m`,
//! additions `dims * m/(w+1) * c_add`, per-dimension tables of `2^(w-2)`
//! points (the `psi` images are charged `c_psi` each instead of a table), and
//! for 4 dimensions the loop `cost(alpha)` once.
//!
//! What this does **not** establish: that a GLV-GLS implementation is faster
//! (no curve, no point arithmetic, no timing); that `lambda` is the action of
//! the loop on points of a constructed `E'` (it is derived from the ring map;
//! the D = -71 Vélu self-test of the earlier round validated the loop-to-
//! eigenvalue correspondence in the prime-field case only).
//!
//! ```bash
//! cargo run --release --example isogeny_loop_glv_gls_lattice -- \
//!     [--bits 127] [--samples 10000] [--seed 20261011] [--max-tries 2000000] \
//!     [--discs -3,-4,-7,-8,-11,-23,-71,-431,-619,-3019] [--json out.json]
//! ```

use crypto_lib::cryptanalysis::curve_construction::weil::is_probable_prime;
use crypto_lib::cryptanalysis::lattice::lll_reduce_hp;
use num_bigint::{BigInt, BigUint, RandBigInt, Sign};
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};
use rand::rngs::StdRng;
use rand::SeedableRng;
use std::fmt::Write as _;

/// Cost per bit of degree in the paper's projective model: min over primes
/// of 7.5 l / log2 l, attained at l = 3.
fn kappa_projective() -> f64 {
    7.5 * 3.0 / 3f64.log2()
}

fn step_cost(l: u64) -> f64 {
    7.5 * l as f64
}

fn kronecker(d: i64, l: u64) -> i32 {
    // Kronecker symbol (d / l) for prime l.
    if l == 2 {
        let r = d.rem_euclid(8);
        return match r {
            1 | 7 => 1,
            3 | 5 => -1,
            _ => 0,
        };
    }
    let dm = d.rem_euclid(l as i64) as u64;
    if dm == 0 {
        return 0;
    }
    let mut res: u64 = 1;
    let mut b = dm;
    let mut e = (l - 1) / 2;
    while e > 0 {
        if e & 1 == 1 {
            res = (res as u128 * b as u128 % l as u128) as u64;
        }
        b = (b as u128 * b as u128 % l as u128) as u64;
        e >>= 1;
    }
    if res == 1 {
        1
    } else {
        -1
    }
}

fn factor_u64(mut n: u64) -> Vec<(u64, u32)> {
    let mut out = vec![];
    let mut q = 2u64;
    while q * q <= n {
        if n.is_multiple_of(q) {
            let mut e = 0;
            while n.is_multiple_of(q) {
                n /= q;
                e += 1;
            }
            out.push((q, e));
        }
        q += if q == 2 { 1 } else { 2 };
    }
    if n > 1 {
        out.push((n, 1));
    }
    out
}

#[derive(Clone, Debug)]
struct Loop {
    a: i64,
    b: i64,
    norm: u64,
    steps: Vec<(u64, u32)>, // (l, number of steps)
    cost: f64,
}

/// d_min of O_K: the smallest norm of an element with B != 0.
fn d_min(d: i64) -> u64 {
    if d.rem_euclid(4) == 0 {
        (-d / 4) as u64
    } else {
        ((1 - d) / 4) as u64
    }
}

/// Every primitive alpha = (A + B sqrt D)/2 in O_K with B >= 1 and
/// N(alpha) <= nmax, as a loop (prime ideal factorization of a primitive
/// element: split l^e gives e steps along one prime, ramified l gives one
/// step; an inert factor cannot divide a primitive element).  Returns the
/// minimum-cost loops (all ties) and the number of elements enumerated.
fn cheapest_loops(d: i64, nmax: u64) -> (Vec<Loop>, u64) {
    let ad = (-d) as u64;
    let mut best: Vec<Loop> = vec![];
    let mut best_cost = f64::INFINITY;
    let mut enumerated = 0u64;
    let even = d.rem_euclid(4) == 0;
    let mut b: i64 = 1;
    while (ad as u128) * (b as u128) * (b as u128) <= 4 * nmax as u128 {
        let amax = ((4 * nmax as u128 - (ad as u128) * (b as u128) * (b as u128)) as f64).sqrt()
            as i64
            + 1;
        for a in -amax..=amax {
            // A must match the parity of B * D so that alpha is in O_K
            if even {
                if a.rem_euclid(2) != 0 {
                    continue;
                }
            } else if (a - b).rem_euclid(2) != 0 {
                continue;
            }
            let n4 = (a as i128) * (a as i128) + (ad as i128) * (b as i128) * (b as i128);
            if n4 % 4 != 0 {
                continue;
            }
            let n = (n4 / 4) as u64;
            if n == 0 || n > nmax {
                continue;
            }
            // primitivity: coordinates in the integral basis {1, omega}
            let (x, y) = if even { (a / 2, b) } else { ((a - b) / 2, b) };
            if x.gcd(&y) != 1 {
                continue;
            }
            enumerated += 1;
            let fac = factor_u64(n);
            let mut ok = true;
            let mut cost = 0.0;
            let mut steps = vec![];
            for &(l, e) in &fac {
                match kronecker(d, l) {
                    -1 => {
                        ok = false;
                        break;
                    }
                    _ => {
                        cost += e as f64 * step_cost(l);
                        steps.push((l, e));
                    }
                }
            }
            if !ok {
                continue;
            }
            if n == 1 {
                cost = 1.0; // a non-trivial automorphism (D = -3, -4): one multiplication
            }
            if cost < best_cost - 1e-9 {
                best_cost = cost;
                best.clear();
            }
            if (cost - best_cost).abs() < 1e-9 {
                best.push(Loop {
                    a,
                    b,
                    norm: n,
                    steps,
                    cost,
                });
            }
        }
        b += 1;
    }
    (best, enumerated)
}

fn modinv(a: &BigInt, m: &BigInt) -> BigInt {
    let e = a.extended_gcd(m);
    assert!(e.gcd.is_one(), "not invertible");
    e.x.mod_floor(m)
}

fn bits(x: &BigInt) -> u64 {
    x.abs().bits()
}

#[derive(Clone, Debug)]
struct Instance {
    p: BigInt,
    t: BigInt,
    v: BigInt,
    h: u64,
    r: BigInt,
    tries: u64,
}

/// Search t, v with p = (t^2 + |D| v^2)/4 prime of exactly `pbits` bits and
/// (p-1)^2 + t^2 = h r with r prime, h <= 64.  Gives up after `max_tries`
/// draws.  For D = -4 no instance exists: with x = t/2 and p = x^2 + v^2,
/// (p-1)^2 + t^2 = |(x+i)^2 + v^2|^2 = (x^2 + (1+v)^2)(x^2 + (1-v)^2), a
/// product of two half-size factors (j = 1728 needs a quartic twist).
fn find_instance(d: i64, pbits: u64, rng: &mut StdRng, max_tries: u64) -> Option<Instance> {
    let ad = BigInt::from(-d);
    let even = d.rem_euclid(4) == 0;
    let four = BigInt::from(4);
    let mut tries = 0u64;
    while tries < max_tries {
        tries += 1;
        // 4p = t^2 + |D| v^2 must reach [2^(pbits+1), 2^(pbits+2)): draw v below
        // sqrt(2^(pbits+2)/|D|) and t below 2^((pbits+2)/2), then reject on p's size.
        let top = BigUint::one() << (pbits + 2);
        let vmax = (&top / BigUint::from((-d) as u64)).sqrt();
        let tmax = top.sqrt();
        let v: BigInt = BigInt::from_biguint(Sign::Plus, rng.gen_biguint_below(&vmax));
        let mut t: BigInt = BigInt::from_biguint(Sign::Plus, rng.gen_biguint_below(&tmax));
        if v.is_zero() {
            continue;
        }
        if even {
            if t.is_odd() {
                t += BigInt::one();
            }
        } else if (&t - &v).is_odd() {
            t += BigInt::one();
        }
        let n4 = &t * &t + &ad * &v * &v;
        if !(&n4 % &four).is_zero() {
            continue;
        }
        let p = &n4 / &four;
        if p.bits() != pbits {
            continue;
        }
        if !is_probable_prime(&p.to_biguint().unwrap()) {
            continue;
        }
        let pm1 = &p - BigInt::one();
        let nt = &pm1 * &pm1 + &t * &t;
        for h in 1u64..=64 {
            let hb = BigInt::from(h);
            if !(&nt % &hb).is_zero() {
                continue;
            }
            let r = &nt / &hb;
            if is_probable_prime(&r.to_biguint().unwrap()) {
                return Some(Instance {
                    p,
                    t,
                    v,
                    h,
                    r,
                    tries,
                });
            }
        }
    }
    None
}

fn det3(m: &[[BigInt; 3]; 3]) -> BigInt {
    &m[0][0] * (&m[1][1] * &m[2][2] - &m[1][2] * &m[2][1])
        - &m[0][1] * (&m[1][0] * &m[2][2] - &m[1][2] * &m[2][0])
        + &m[0][2] * (&m[1][0] * &m[2][1] - &m[1][1] * &m[2][0])
}

/// Babai rounding data for the target (k, 0, ..., 0): c_j = round(k * C_j0 / det)
/// where C_j0 is the cofactor of entry (j, 0) of the row basis.
struct Babai {
    basis: Vec<Vec<BigInt>>,
    cof: Vec<BigInt>,
    det: BigInt,
}

fn babai_setup(basis: Vec<Vec<BigInt>>) -> Babai {
    let n = basis.len();
    let (cof, det) = if n == 2 {
        // C_00 = b11, C_10 = -b01
        let cof = vec![basis[1][1].clone(), -basis[0][1].clone()];
        let det = &basis[0][0] * &basis[1][1] - &basis[0][1] * &basis[1][0];
        (cof, det)
    } else {
        assert_eq!(n, 4);
        let mut cof = vec![];
        for j in 0..4 {
            let rows: Vec<usize> = (0..4).filter(|&i| i != j).collect();
            let mut m: [[BigInt; 3]; 3] = Default::default();
            for (ri, &i) in rows.iter().enumerate() {
                for c in 1..4 {
                    m[ri][c - 1] = basis[i][c].clone();
                }
            }
            let s = if j % 2 == 0 {
                BigInt::one()
            } else {
                -BigInt::one()
            };
            cof.push(s * det3(&m));
        }
        let det: BigInt = (0..4).map(|j| &basis[j][0] * &cof[j]).sum();
        (cof, det)
    };
    assert!(!det.is_zero());
    Babai { basis, cof, det }
}

fn round_div(num: &BigInt, den: &BigInt) -> BigInt {
    let (num, den) = if den.is_negative() {
        (-num, -den)
    } else {
        (num.clone(), den.clone())
    };
    let two = BigInt::from(2);
    (&two * num + &den).div_floor(&(&two * den))
}

fn decompose(bb: &Babai, k: &BigInt) -> Vec<BigInt> {
    let n = bb.basis.len();
    let mut out = vec![BigInt::zero(); n];
    out[0] = k.clone();
    for j in 0..n {
        let c = round_div(&(k * &bb.cof[j]), &bb.det);
        for i in 0..n {
            out[i] -= &c * &bb.basis[j][i];
        }
    }
    out
}

#[derive(Default)]
struct DecompStats {
    max_bits: u64,
    mean_bits: f64,
    basis_max_bits: u64,
}

fn measure(
    basis: Vec<Vec<BigInt>>,
    eig: &[BigInt],
    r: &BigInt,
    samples: u64,
    rng: &mut StdRng,
) -> DecompStats {
    let mut b = basis;
    lll_reduce_hp(&mut b, 0.99).expect("LLL");
    let basis_max_bits = b.iter().flat_map(|row| row.iter().map(bits)).max().unwrap();
    // every reduced vector must still lie in the lattice
    for row in &b {
        let s: BigInt = row.iter().zip(eig).map(|(x, e)| x * e).sum();
        assert!(s.mod_floor(r).is_zero(), "LLL output left the lattice");
    }
    let bb = babai_setup(b);
    let mut max_bits = 0u64;
    let mut sum_bits = 0f64;
    let ru = r.to_biguint().unwrap();
    for _ in 0..samples {
        let k = BigInt::from_biguint(Sign::Plus, rng.gen_biguint_below(&ru));
        let dvec = decompose(&bb, &k);
        let s: BigInt = dvec.iter().zip(eig).map(|(x, e)| x * e).sum();
        assert!(
            (s - &k).mod_floor(r).is_zero(),
            "decomposition congruence failed"
        );
        let m = dvec.iter().map(bits).max().unwrap();
        max_bits = max_bits.max(m);
        sum_bits += m as f64;
    }
    DecompStats {
        max_bits,
        mean_bits: sum_bits / samples as f64,
        basis_max_bits,
    }
}

struct Model {
    w: u32,
    c_dbl: f64,
    c_add: f64,
    c_psi: f64,
}

impl Model {
    fn table(&self) -> f64 {
        (1u64 << (self.w - 2)) as f64 * self.c_add
    }
    fn gls2(&self, m: f64) -> f64 {
        // table for P, psi applied to the table
        self.c_dbl * m
            + 2.0 * m / (self.w as f64 + 1.0) * self.c_add
            + self.table()
            + (1u64 << (self.w - 2)) as f64 * self.c_psi
    }
    fn glvgls4(&self, m: f64, loop_cost: f64) -> f64 {
        // tables for P and alpha(P), psi applied to both tables
        self.c_dbl * m
            + 4.0 * m / (self.w as f64 + 1.0) * self.c_add
            + 2.0 * self.table()
            + 2.0 * (1u64 << (self.w - 2)) as f64 * self.c_psi
            + loop_cost
    }
}

fn arg<T: std::str::FromStr>(args: &[String], name: &str, default: T) -> T {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let pbits: u64 = arg(&args, "--bits", 127);
    let samples: u64 = arg(&args, "--samples", 10_000);
    let seed: u64 = arg(&args, "--seed", 20_261_011);
    let max_tries: u64 = arg(&args, "--max-tries", 2_000_000);
    let discs: Vec<i64> = arg(
        &args,
        "--discs",
        "-3,-4,-7,-8,-11,-23,-71,-431,-619,-3019".to_string(),
    )
    .split(',')
    .map(|s| s.trim().parse().expect("disc"))
    .collect();
    let json_out: Option<String> = args
        .iter()
        .position(|a| a == "--json")
        .and_then(|i| args.get(i + 1))
        .cloned();
    let model = Model {
        w: 4,
        c_dbl: 8.0,
        c_add: 12.0,
        c_psi: 4.0,
    };
    let kappa = kappa_projective();

    println!(
        "# Isogeny-loop endomorphisms in 4-dimensional GLV-GLS (decomposition-lattice level)\n"
    );
    println!(
        "p bits {pbits}, samples per lattice {samples}, seed {seed}; projective cost model 7.5 l M per l-step, \
         kappa = {kappa:.4} M per bit of degree; operation-count model w = {}, c_dbl = {}, c_add = {}, c_psi = {} (F_(p^2) multiplications)\n",
        model.w, model.c_dbl, model.c_add, model.c_psi
    );

    let mut table = String::new();
    let mut loops_table = String::new();
    let mut json_rows = vec![];
    writeln!(loops_table, "| D | d_min | cheapest loop alpha | N(alpha) | steps (l: count) | cost M | kappa log2 d_min (bound) | elements enumerated | norm bound |").unwrap();
    writeln!(loops_table, "|--:|--:|:--|--:|:--|--:|--:|--:|--:|").unwrap();
    writeln!(
        table,
        "| D | p bits | twist cofactor h | log2 r | tries | embedding | GLS-2 max/mean bits (log2 r/2) | GLV-GLS-4 max/mean bits (log2 r/4) | GLV-GLS-4 basis max bits | model GLS-2 M | model GLV-GLS-4 M | model ratio GLS-2/GLV-GLS-4 |"
    )
    .unwrap();
    writeln!(table, "|--:|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|").unwrap();

    for (di, &d) in discs.iter().enumerate() {
        let mut rng = StdRng::seed_from_u64(seed ^ ((di as u64 + 1) * 0x9E37_79B9_7F4A_7C15));
        // 1. certified cheapest loop: first pass finds an upper bound, second pass certifies
        let dm = d_min(d);
        let (first, _) = cheapest_loops(d, dm.max(1) * 64);
        let ub = first.first().map(|l| l.cost).unwrap_or(f64::INFINITY);
        let nmax = if ub.is_finite() {
            (2f64.powf(ub / kappa)).floor() as u64
        } else {
            dm * 4096
        };
        let nmax = nmax.max(dm.max(1) * 64);
        let (best, enumerated) = cheapest_loops(d, nmax);
        let lp = best.first().expect("a loop exists").clone();
        let bound = if dm > 1 {
            kappa * (dm as f64).log2()
        } else {
            0.0
        };
        assert!(
            lp.cost + 1e-9 >= bound || lp.norm == 1,
            "cost below the degree bound"
        );
        let steps_s = lp
            .steps
            .iter()
            .map(|(l, e)| format!("{l}: {e}"))
            .collect::<Vec<_>>()
            .join(", ");
        writeln!(
            loops_table,
            "| {d} | {dm} | ({} + {} sqrt({d}))/2 (+{} ties) | {} | {} | {:.1} | {:.1} | {} | {} |",
            lp.a,
            lp.b,
            best.len() - 1,
            lp.norm,
            if steps_s.is_empty() {
                "automorphism".into()
            } else {
                steps_s
            },
            lp.cost,
            bound,
            enumerated,
            nmax
        )
        .unwrap();

        // 2. GLS instance
        let Some(inst) = find_instance(d, pbits, &mut rng, max_tries) else {
            writeln!(
                table,
                "| {d} | {pbits} | — | — | {max_tries} | — | no instance within {max_tries} tries | — | — | — | — | — |"
            )
            .unwrap();
            json_rows.push(serde_json::json!({"D": d, "d_min": dm, "instance": null, "max_tries": max_tries,
                "loop": {"A": lp.a, "B": lp.b, "norm": lp.norm, "steps": lp.steps, "cost_M": lp.cost}}));
            continue;
        };
        let (p, t, v, r) = (&inst.p, &inst.t, &inst.v, &inst.r);
        let pm1 = p - BigInt::one();
        let mu = (&pm1 * modinv(t, r)).mod_floor(r);
        assert!(
            (&mu * &mu + BigInt::one()).mod_floor(r).is_zero(),
            "mu^2 != -1"
        );
        assert!(
            (&mu * &mu - t * &mu + p).mod_floor(r).is_zero(),
            "mu not a Frobenius eigenvalue"
        );
        let s = ((BigInt::from(2) * &mu - t) * modinv(v, r)).mod_floor(r);
        assert!(
            (&s * &s - BigInt::from(d)).mod_floor(r).is_zero(),
            "s^2 != D"
        );
        let inv2 = modinv(&BigInt::from(2), r);
        let a_big = BigInt::from(lp.a);
        let b_big = BigInt::from(lp.b);
        let n_big = BigInt::from(lp.norm);
        let lr = (r.bits()) as f64;

        // GLS-2 lattice (independent of alpha)
        let gls2 = measure(
            vec![
                vec![r.clone(), BigInt::zero()],
                vec![(-&mu).mod_floor(r), BigInt::one()],
            ],
            &[BigInt::one(), mu.clone()],
            r,
            samples,
            &mut rng,
        );

        for (emb, sgn) in [("s", 1i32), ("-s", -1i32)] {
            let ss = if sgn > 0 {
                s.clone()
            } else {
                (-&s).mod_floor(r)
            };
            let lam = ((&a_big + &b_big * &ss) * &inv2).mod_floor(r);
            assert!(
                (&lam * &lam - &a_big * &lam + &n_big)
                    .mod_floor(r)
                    .is_zero(),
                "lambda not a root of the char poly"
            );
            let lm = (&lam * &mu).mod_floor(r);
            let eig = [BigInt::one(), mu.clone(), lam.clone(), lm.clone()];
            let basis = vec![
                vec![r.clone(), BigInt::zero(), BigInt::zero(), BigInt::zero()],
                vec![
                    (-&mu).mod_floor(r),
                    BigInt::one(),
                    BigInt::zero(),
                    BigInt::zero(),
                ],
                vec![
                    (-&lam).mod_floor(r),
                    BigInt::zero(),
                    BigInt::one(),
                    BigInt::zero(),
                ],
                vec![
                    (-&lm).mod_floor(r),
                    BigInt::zero(),
                    BigInt::zero(),
                    BigInt::one(),
                ],
            ];
            let g4 = measure(basis, &eig, r, samples, &mut rng);
            let c2 = model.gls2(gls2.mean_bits);
            let c4 = model.glvgls4(g4.mean_bits, lp.cost);
            writeln!(
                table,
                "| {d} | {pbits} | {} | {:.1} | {} | {emb} | {}/{:.1} ({:.1}) | {}/{:.1} ({:.1}) | {} | {:.0} | {:.0} | {:.3} |",
                inst.h,
                lr,
                inst.tries,
                gls2.max_bits,
                gls2.mean_bits,
                lr / 2.0,
                g4.max_bits,
                g4.mean_bits,
                lr / 4.0,
                g4.basis_max_bits,
                c2,
                c4,
                c2 / c4
            )
            .unwrap();
            json_rows.push(serde_json::json!({
                "D": d, "d_min": dm,
                "loop": {"A": lp.a, "B": lp.b, "norm": lp.norm, "steps": lp.steps, "cost_M": lp.cost,
                          "ties": best.len() - 1, "enumerated": enumerated, "norm_bound": nmax,
                          "degree_bound_M": bound},
                "instance": {"p": p.to_string(), "t": t.to_string(), "v": v.to_string(),
                              "h": inst.h, "r": r.to_string(), "log2_r": lr, "tries": inst.tries},
                "mu": mu.to_string(), "embedding": emb, "lambda": lam.to_string(),
                "gls2": {"max_bits": gls2.max_bits, "mean_bits": gls2.mean_bits, "basis_max_bits": gls2.basis_max_bits},
                "glvgls4": {"max_bits": g4.max_bits, "mean_bits": g4.mean_bits, "basis_max_bits": g4.basis_max_bits},
                "model": {"w": model.w, "c_dbl": model.c_dbl, "c_add": model.c_add, "c_psi": model.c_psi,
                           "gls2_M": c2, "glvgls4_M": c4, "ratio": c2 / c4},
            }));
        }
    }

    println!("## Cheapest loop per discriminant (exhaustive up to the norm bound)\n");
    print!("{loops_table}");
    println!("\n## Decomposition lattices on the GLS twist E'(F_(p^2))\n");
    print!("{table}");
    println!(
        "\nAll identities asserted per instance: mu^2 = -1, mu^2 - t mu + p = 0, s^2 = D, lambda^2 - Tr(alpha) lambda + N(alpha) = 0 (mod r); \
         every LLL vector in the lattice; every decomposition congruent to its scalar. Model columns are derived operation counts, not measurements."
    );
    if let Some(path) = json_out {
        let doc = serde_json::json!({
            "schema": "isogeny_loop_glv_gls_lattice/v1",
            "command": args.join(" "),
            "bits": pbits, "samples": samples, "seed": seed, "kappa": kappa,
            "rows": json_rows,
        });
        std::fs::write(&path, serde_json::to_string_pretty(&doc).unwrap()).expect("write json");
        eprintln!("wrote {path}");
    }
    let _ = BigUint::zero().to_u64();
}
