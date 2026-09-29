//! End-to-end outputs of the cryptanalysis walks that moved to the
//! variable-time point arithmetic (`Point::add_vartime` and siblings),
//! pinned to what the constant-time arithmetic produced before the move.
//!
//! The unit tests in `ecc::point` compare the vartime operations with
//! the constant-time ones one call at a time.  These run whole solves —
//! every `ecdlp_variants` solver, the collaborative rho walkers with and
//! without the negation map, Pohlig–Hellman and Cheon — over several
//! curves, targets, seeds and parameters, and hash everything they
//! report: the answer, the charged counters (`group_ops`,
//! `field_inversions`, `table_size`, baby steps, exponentiations) and
//! the distinguished points.  Each digest was produced by the tree at
//! b072fcf5, before the change; a mismatch means a walk now takes a
//! different path or charges a different count.
//!
//! The second half holds the vartime operations to the constant-time
//! ones off the prime-field domain: every residue of small moduli, and
//! random operations on composite moduli and unreduced coordinates.

use crypto_lib::cryptanalysis::cheon_attack::cheon_attack;
use crypto_lib::cryptanalysis::ecdlp_variants::{
    bsgs_average_case, bsgs_interleaving, bsgs_interleaving_block, bsgs_interleaving_negation,
    bsgs_negation, bsgs_textbook, demo_group_mid, demo_group_small, gaudry_schost,
    gaudry_schost_montgomery, gaudry_schost_negation, grumpy_giants, grumpy_giants_block,
    grumpy_giants_negation, DlpSolution, EcGroup, GaudrySchostOptions,
};
use crypto_lib::cryptanalysis::pohlig_hellman::pohlig_hellman_curve;
use crypto_lib::cryptanalysis::pollard_collab::{demo_curve, run_walker, JobSpec, WalkerOutcome};
use crypto_lib::ecc::{CurveParams, Point};
use num_bigint::BigUint;
use num_traits::Zero;
use std::fmt::Write;

/// FNV-1a over a transcript: enough to pin it, short enough to embed.
fn fnv(s: &str) -> u64 {
    s.bytes().fold(0xcbf2_9ce4_8422_2325, |h, b| {
        (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3)
    })
}

fn sol(out: &mut String, tag: &str, r: Result<DlpSolution, &'static str>) {
    match r {
        Ok(s) => writeln!(
            out,
            "{tag} ok {} {} {} {}",
            s.x, s.group_ops, s.field_inversions, s.table_size
        ),
        Err(e) => writeln!(out, "{tag} err {e}"),
    }
    .unwrap();
}

fn pt(p: &Point) -> String {
    match p {
        Point::Infinity => "O".into(),
        Point::Affine { x, y } => format!("({},{})", x.value, y.value),
    }
}

fn ecdlp_transcript(name: &str, group: &EcGroup) -> String {
    let mut out = String::new();
    let n = group.order().clone();
    let one = BigUint::from(1u32);
    let mut ks = vec![one.clone(), BigUint::from(2u32), &n - &one, &n >> 1];
    let mut s = 0x243f_6a88_85a3_08d3u64;
    for _ in 0..3 {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        ks.push(BigUint::from(s) % &n);
    }
    ks.push(BigUint::zero());
    for k in &ks {
        let q = group.mul_setup(k);
        let t = format!("{name} k={k}");
        sol(&mut out, &format!("{t} textbook"), bsgs_textbook(group, &q));
        sol(
            &mut out,
            &format!("{t} average"),
            bsgs_average_case(group, &q),
        );
        sol(
            &mut out,
            &format!("{t} interleave"),
            bsgs_interleaving(group, &q),
        );
        sol(&mut out, &format!("{t} grumpy"), grumpy_giants(group, &q));
        sol(&mut out, &format!("{t} neg"), bsgs_negation(group, &q));
        sol(
            &mut out,
            &format!("{t} interleave-neg"),
            bsgs_interleaving_negation(group, &q),
        );
        sol(
            &mut out,
            &format!("{t} grumpy-neg"),
            grumpy_giants_negation(group, &q),
        );
        for b in [1usize, 7, 32] {
            sol(
                &mut out,
                &format!("{t} interleave-block{b}"),
                bsgs_interleaving_block(group, &q, b),
            );
            sol(
                &mut out,
                &format!("{t} grumpy-block{b}"),
                grumpy_giants_block(group, &q, b),
            );
        }
        for (seed, dp_bits, jumps, block) in
            [(1u64, 2u8, 16usize, 1usize), (2, 4, 32, 8), (3, 6, 7, 32)]
        {
            let opts = GaudrySchostOptions {
                dp_bits,
                num_jumps: jumps,
                max_walkers: 1 << 16,
                max_steps_per_walker: 0,
                block,
                seed: Some(seed),
            };
            let o = format!("{t} gs seed={seed} dp={dp_bits}");
            sol(
                &mut out,
                &format!("{o} plain"),
                gaudry_schost(group, &q, &opts),
            );
            sol(
                &mut out,
                &format!("{o} neg"),
                gaudry_schost_negation(group, &q, &opts),
            );
            sol(
                &mut out,
                &format!("{o} mont"),
                gaudry_schost_montgomery(group, &q, &opts),
            );
        }
    }
    out
}

#[test]
fn ecdlp_variants_solves_are_unchanged() {
    let small = ecdlp_transcript("small", &demo_group_small());
    let mid = ecdlp_transcript("mid", &demo_group_mid());
    assert_eq!(
        (fnv(&small), fnv(&mid)),
        (ECDLP_SMALL, ECDLP_MID),
        "{small}\n{mid}"
    );
}

fn collab_transcript(curve_name: &str, walkers: u64) -> String {
    let mut out = String::new();
    let curve = demo_curve(curve_name).expect("demo curve");
    let q = curve
        .generator()
        .scalar_mul(&(&curve.n * 2u32 / 3u32), &curve.a_fe());
    for negation in [false, true] {
        for extra_dp in [0u8, 2] {
            let mut spec = JobSpec::new(&curve, &q, "pinned", 0x5eed_0000 + u64::from(extra_dp))
                .expect("job spec");
            spec.negation_map = negation;
            spec.dp_bits += extra_dp;
            let ctx = spec.build().expect("job context");
            write!(out, "{curve_name} neg={negation} dp={}:", spec.dp_bits).unwrap();
            for b in &ctx.branches {
                write!(out, " {}", pt(&b.point)).unwrap();
            }
            writeln!(out).unwrap();
            for i in 0..walkers {
                match run_walker(&ctx, i) {
                    WalkerOutcome::Dp(r) => writeln!(
                        out,
                        "  {i} dp {} {} {} {} {} {}",
                        r.walker, r.steps, r.x, r.y, r.a, r.b
                    ),
                    WalkerOutcome::DeadTrail { steps } => writeln!(out, "  {i} dead {steps}"),
                }
                .unwrap();
            }
        }
    }
    out
}

#[test]
fn collab_walkers_are_unchanged() {
    let t = [
        collab_transcript("demo-small", 24),
        collab_transcript("demo-mid", 24),
        collab_transcript("demo-32", 6),
    ];
    let got = [fnv(&t[0]), fnv(&t[1]), fnv(&t[2])];
    assert_eq!(got, COLLAB, "{}", t.concat());
}

/// `y² = x³ + 4` over a 47-bit prime; order `3 · 19⁴ · 151 · 181 · 13171`
/// (the perfbench Pohlig–Hellman curve).
fn ph_curve() -> CurveParams {
    CurveParams {
        name: "pinned-j0-smooth47",
        p: BigUint::from(140_737_529_461_027u64),
        a: BigUint::from(0u32),
        b: BigUint::from(4u32),
        gx: BigUint::from(16_535_943_672_630u64),
        gy: BigUint::from(38_148_592_391_465u64),
        n: BigUint::from(140_737_531_856_763u64),
        h: 1,
    }
}

#[test]
fn pohlig_hellman_and_cheon_are_unchanged() {
    let mut out = String::new();
    let curves = [
        ph_curve(),
        demo_curve("demo-small").unwrap(),
        demo_curve("demo-mid").unwrap(),
    ];
    for curve in &curves {
        let g = curve.generator();
        let a = curve.a_fe();
        for d in [0u64, 1, 2, 12_345, 99_991, 140_737_531_856_000] {
            let d = BigUint::from(d) % &curve.n;
            let q = g.scalar_mul(&d, &a);
            for bound in [200u64, 1 << 14, 1 << 17] {
                let r = pohlig_hellman_curve(curve, &g, &q, &curve.n, bound);
                writeln!(
                    out,
                    "{} d={d} bound={bound}: {:?} {:?} {:?} {}",
                    curve.name, r.recovered_d, r.factors, r.residues, r.total_steps
                )
                .unwrap();
            }
            // A mis-stated order: the searches run on the wrong subgroups.
            let wrong = &curve.n * 2u32;
            let r = pohlig_hellman_curve(curve, &g, &q, &wrong, 1 << 14);
            writeln!(
                out,
                "  wrong-order: {:?} {:?} {}",
                r.recovered_d, r.residues, r.total_steps
            )
            .unwrap();
        }
    }
    // Cheon on its toy curve of order 199 (199 − 1 = 2 · 3² · 11).
    let curve = CurveParams {
        name: "cheon-199",
        p: BigUint::from(211u32),
        a: BigUint::zero(),
        b: BigUint::from(2u32),
        gx: BigUint::from(4u32),
        gy: BigUint::from(53u32),
        n: BigUint::from(199u32),
        h: 1,
    };
    let g = curve.generator();
    let a = curve.a_fe();
    for d in [1u32, 2, 3, 6, 9, 11, 18, 22, 33, 66, 99, 198, 7] {
        for alpha in [2u32, 5, 73, 100, 197] {
            let alpha_b = BigUint::from(alpha);
            let ag = g.scalar_mul(&alpha_b, &a);
            let adg = g.scalar_mul(&alpha_b.modpow(&BigUint::from(d), &curve.n), &a);
            let r = cheon_attack(&curve, &g, &ag, &adg, &curve.n, &BigUint::from(d));
            writeln!(
                out,
                "cheon d={d} alpha={alpha}: {:?} {} {} {:?}",
                r.recovered_alpha, r.step1_exps, r.step2_exps, r.error
            )
            .unwrap();
        }
    }
    assert_eq!(fnv(&out), PH_CHEON, "{out}");
}

const ECDLP_SMALL: u64 = 12_107_984_968_993_051_905;
const ECDLP_MID: u64 = 3_865_584_706_469_309_340;
const COLLAB: [u64; 3] = [
    13_772_181_218_826_505_772,
    6_695_007_987_417_753_287,
    9_866_830_489_845_429_253,
];
const PH_CHEON: u64 = 16_994_138_649_881_210_032;

// ── Off the prime-field domain ───────────────────────────────────────
//
// `FieldElement` is documented for a prime modulus, but a malformed job
// file can hand the walks any modulus, and the fields are public, so a
// caller can build an unreduced coordinate.  Over a prime the vartime
// operations must agree with the constant-time ones exactly, on any
// input.  Over a composite modulus they may return a different point (a
// true inverse where Fermat's `a^(p−2)` is none), but on reduced
// coordinates — all a job file can produce, since it goes through
// `FieldElement::new` — they must panic on exactly the inputs the
// constant-time ones panic on.  (A composite modulus *and* unreduced
// coordinates together is not covered: there the panics are those of
// a `BigUint` underflow on values that already differ.)

use crypto_lib::ecc::FieldElement;
use num_integer::Integer;
use std::panic::{catch_unwind, AssertUnwindSafe};

fn elem(v: BigUint, m: &BigUint) -> FieldElement {
    FieldElement {
        value: v,
        modulus: m.clone(),
    }
}

/// `f()`'s point with both moduli, or `None` if it panicked.
fn run(f: impl FnOnce() -> Point) -> Option<String> {
    catch_unwind(AssertUnwindSafe(f)).ok().map(|p| match p {
        Point::Infinity => "O".into(),
        Point::Affine { x, y } => format!("({}/{},{}/{})", x.value, x.modulus, y.value, y.modulus),
    })
}

/// Every residue of every modulus up to 700, and some unreduced ones,
/// then full-word moduli: a unit gets its inverse, which is `inv`'s
/// whenever the modulus is prime; a non-unit gets `inv`'s value; zero
/// gets `None` from both.
#[test]
fn inv_vartime_on_every_small_modulus() {
    for m in 2u64..=700 {
        let mb = BigUint::from(m);
        let prime = (2..m).take_while(|d| d * d <= m).all(|d| m % d != 0);
        for a in 0..(m + m / 3) {
            let x = elem(BigUint::from(a), &mb);
            let (s, t) = (x.inv_vartime(), x.inv());
            if a == 0 {
                assert!(s.is_none() && t.is_none());
                continue;
            }
            let (s, t) = (s.expect("nonzero"), t.expect("nonzero"));
            assert_eq!(s.modulus, mb);
            if a.gcd(&m) == 1 {
                assert!(s.value < mb, "{a}^-1 mod {m}");
                assert_eq!((&s.value * a) % m, BigUint::from(1u32), "{a}^-1 mod {m}");
                if prime {
                    assert_eq!(s.value, t.value, "{a}^-1 mod {m}");
                }
            } else {
                assert_eq!(s.value, t.value, "non-unit {a} mod {m}");
            }
        }
    }
    let mut r = 0x0123_4567_89ab_cdefu64;
    for m in [
        u64::MAX,
        u64::MAX - 1,
        1u64 << 63,
        (1u64 << 63) + 1,
        u64::MAX - 58,
        4_294_967_295,
        4_294_967_291 * 3,
    ] {
        let mb = BigUint::from(m);
        for i in 0..300u64 {
            r ^= r << 13;
            r ^= r >> 7;
            r ^= r << 17;
            let a = match i {
                0 => 1,
                1 => m - 1,
                _ => r % m,
            };
            if a == 0 {
                continue;
            }
            let x = elem(BigUint::from(a), &mb);
            let s = x.inv_vartime().unwrap();
            if a.gcd(&m) == 1 {
                assert_eq!((&s.value * a) % &mb, BigUint::from(1u32), "{a}^-1 mod {m}");
            } else {
                assert_eq!(s.value, x.inv().unwrap().value, "non-unit {a} mod {m}");
            }
        }
    }
}

/// Random `add`, `double`, `neg` and `scalar_mul` on prime and composite
/// moduli, one word and several, with the special cases drawn on
/// purpose: over a prime the vartime and constant-time operations
/// return the same point or both panic, on reduced and unreduced
/// coordinates alike; over a composite modulus, on reduced
/// coordinates, they panic on the same inputs.
#[test]
fn vartime_ops_panic_where_constant_time_ops_do() {
    use num_bigint::RandBigInt;
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};
    let big = |s: &str| BigUint::parse_bytes(s.as_bytes(), 10).unwrap();
    let moduli: Vec<(BigUint, bool)> = vec![
        (big("2"), true),
        (big("3"), true),
        (big("97"), true),
        (big("65521"), true),
        (big("18446744073709551557"), true),
        (big("340282366920938463463374607431768211297"), true),
        (big("4"), false),
        (big("9"), false),
        (big("15"), false),
        (big("91"), false),
        (big("2021027"), false),
        (big("4294967295"), false),
        (big("18446744073709551615"), false),
        (big("18446744073709551614"), false),
        (big("340282366920938463463374607431768211456"), false),
        (
            &big("2147483647") * &big("618970019642690137449562111"),
            false,
        ),
    ];
    let hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let mut rng = StdRng::seed_from_u64(0x0dd_ba11);
    let mut bad = Vec::new();
    for (m, prime) in &moduli {
        // Unreduced coordinates (`v + m`, `v + 2m`) only where the
        // modulus is prime.
        let lift = if *prime { 16 } else { 0 };
        for _ in 0..1000 {
            let coord = |rng: &mut StdRng| -> BigUint {
                let v = match rng.gen_range(0..10) {
                    0 => BigUint::zero(),
                    1 => BigUint::from(1u32),
                    2 => m - 1u32,
                    _ => rng.gen_biguint_below(m),
                };
                match rng.gen_range(0..16 + lift) {
                    0 if lift > 0 => v + m,
                    1 if lift > 0 => v + m + m,
                    _ => v,
                }
            };
            let a = elem(coord(&mut rng), m);
            let p1 = if rng.gen_range(0..12) == 0 {
                Point::Infinity
            } else {
                Point::Affine {
                    x: elem(coord(&mut rng), m),
                    y: elem(coord(&mut rng), m),
                }
            };
            // The same point, its negative, the same `x` with another
            // `y`, or unrelated.
            let p2 = match (rng.gen_range(0..5), &p1) {
                (0, _) => p1.clone(),
                // (`−y` from the residue: `neg` itself underflows on an
                // unreduced `y`.)
                (1, Point::Affine { x, y }) => Point::Affine {
                    x: x.clone(),
                    y: elem((m - &y.value % m) % m, m),
                },
                (2, Point::Affine { x, .. }) => Point::Affine {
                    x: x.clone(),
                    y: elem(coord(&mut rng), m),
                },
                _ => Point::Affine {
                    x: elem(coord(&mut rng), m),
                    y: elem(coord(&mut rng), m),
                },
            };
            let k = match rng.gen_range(0..3) {
                0 => BigUint::from(rng.gen_range(0..5u32)),
                1 => BigUint::from(rng.gen::<u64>()),
                _ => {
                    let bits = rng.gen_range(1..150);
                    rng.gen_biguint(bits)
                }
            };
            let pairs = [
                (
                    "add",
                    run(|| p1.add_vartime(&p2, &a)),
                    run(|| p1.add(&p2, &a)),
                ),
                ("dbl", run(|| p1.double_vartime(&a)), run(|| p1.double(&a))),
                ("neg", run(|| p1.neg_vartime()), run(|| p1.neg())),
                (
                    "mul",
                    run(|| p1.scalar_mul_vartime(&k, &a)),
                    run(|| p1.scalar_mul(&k, &a)),
                ),
            ];
            for (op, s, t) in pairs {
                if s.is_some() != t.is_some() || (*prime && s != t) {
                    bad.push(format!(
                        "{op} mod {m}: {p1:?} {p2:?} a={a} k={k}: vartime {s:?}, constant-time {t:?}"
                    ));
                }
            }
        }
    }
    std::panic::set_hook(hook);
    assert!(
        bad.is_empty(),
        "{} mismatches, first:\n{}",
        bad.len(),
        bad[0]
    );
}
