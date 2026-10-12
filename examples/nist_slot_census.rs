//! NIST prime-curve trapdoor-slot census.
//!
//! Executes, for P-192, P-224, P-256, P-384 and P-521, every check in
//! `research/nist_trapdoor_slot_program_20261011/PROTOCOL.md` that is
//! computable from the public parameters alone (slots S2, S3, S4, S5, S6
//! of that protocol).  Nothing here solves a discrete logarithm; the
//! output is a table of invariants and of the degrees of freedom the
//! X9.62 generation procedure leaves open.
//!
//! ```text
//! cargo run --release --example nist_slot_census            # Markdown
//! cargo run --release --example nist_slot_census -- --json  # machine form
//! ```

use crypto_lib::cryptanalysis::p256_backdoor_map::{verify_seed_derivation, x962_candidate_r};
use crypto_lib::cryptanalysis::p256_structural::{
    curve_structural_report, deep_factor, trial_factor,
};
use crypto_lib::ecc::curve::CurveParams;
use crypto_lib::hash::sha1::sha1;
use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::{One, ToPrimitive, Zero};
use serde::Serialize;

/// Published FIPS 186-4 / SEC 2 seeds (Appendix D.1.2 / SEC 2 §2.4).
fn seed_for(name: &str) -> [u8; 20] {
    let hex_seed = match name {
        "P-192" => "3045ae6fc8422f64ed579528d38120eae12196d5",
        "P-224" => "bd71344799d5c7fcdc45b59fa3b9ab8f6a948bc5",
        "P-256" => "c49d360886e704936a6678e1139d26b7819f7e90",
        "P-384" => "a335926aa319a27a1d00896a6773a4827acdac73",
        "P-521" => "d09e8800291cb85396cc6717393284aaa0da64ba",
        other => panic!("no published seed for {other}"),
    };
    let bytes = hex::decode(hex_seed).expect("seed hex");
    let mut out = [0u8; 20];
    out.copy_from_slice(&bytes);
    out
}

fn legendre(x: &BigUint, p: &BigUint) -> i8 {
    if x.is_zero() {
        return 0;
    }
    let e = (p - BigUint::one()) >> 1u32;
    let r = x.modpow(&e, p);
    if r.is_one() {
        1
    } else {
        -1
    }
}

/// Tonelli–Shanks.  Returns the root with the smaller integer value.
fn sqrt_mod(x: &BigUint, p: &BigUint) -> Option<BigUint> {
    if x.is_zero() {
        return Some(BigUint::zero());
    }
    if legendre(x, p) != 1 {
        return None;
    }
    let one = BigUint::one();
    let two = BigUint::from(2u32);
    let pm1 = p - &one;
    let root = if (p % 4u32) == BigUint::from(3u32) {
        x.modpow(&((p + &one) >> 2u32), p)
    } else {
        let mut q = pm1.clone();
        let mut s = 0u32;
        while q.is_even() {
            q >>= 1u32;
            s += 1;
        }
        let mut z = two.clone();
        while legendre(&z, p) != -1 {
            z += &one;
        }
        let mut m = s;
        let mut c = z.modpow(&q, p);
        let mut t = x.modpow(&q, p);
        let mut r = x.modpow(&((&q + &one) >> 1u32), p);
        while !t.is_one() {
            let mut i = 0u32;
            let mut tt = t.clone();
            while !tt.is_one() {
                tt = (&tt * &tt) % p;
                i += 1;
            }
            let b = c.modpow(&two.pow(m - i - 1), p);
            m = i;
            c = (&b * &b) % p;
            t = (&t * &c) % p;
            r = (&r * &b) % p;
        }
        r
    };
    let other = p - &root;
    Some(if root < other { root } else { other })
}

fn mod_inverse(x: &BigUint, p: &BigUint) -> BigUint {
    x.modpow(&(p - BigUint::from(2u32)), p)
}

/// Multiplicative order of `p` modulo `n` if it is at most `max_k`.
fn embedding_degree_bound(p: &BigUint, n: &BigUint, max_k: u32) -> Option<u32> {
    let base = p % n;
    let mut acc = BigUint::one();
    for k in 1..=max_k {
        acc = (&acc * &base) % n;
        if acc.is_one() {
            return Some(k);
        }
    }
    None
}

/// Largest square factor found by trial division to `bound`, and the
/// cofactor's status: the residue `|D| / smooth` is reported as square,
/// probable prime, or composite-unknown.
#[derive(Serialize, Clone)]
struct DiscriminantProfile {
    abs_disc_bits: u64,
    trial_bound: u64,
    small_factors: Vec<(u64, u32)>,
    square_part_bits_from_small_factors: u64,
    residue_bits: u64,
    residue_is_perfect_square: bool,
    residue_is_probable_prime: bool,
}

fn discriminant_profile(d: &BigUint, bound: u64) -> DiscriminantProfile {
    let (factors, residue) = trial_factor(d, bound);
    let mut square_part = BigUint::one();
    for (q, e) in &factors {
        if *e >= 2 {
            square_part *= BigUint::from(*q).pow(e - (e % 2));
        }
    }
    let root_residue = residue.sqrt();
    let residue_is_perfect_square = &root_residue * &root_residue == residue;
    DiscriminantProfile {
        abs_disc_bits: d.bits(),
        trial_bound: bound,
        small_factors: factors,
        square_part_bits_from_small_factors: square_part.bits(),
        residue_bits: residue.bits(),
        residue_is_perfect_square,
        residue_is_probable_prime: residue > BigUint::one()
            && crypto_lib::asymmetric::rsa::is_prime(&residue),
    }
}

/// Slot S3: the X9.62 root-sign freedom.
#[derive(Serialize, Clone)]
struct RootSignFreedom {
    seed_derivation_holds: bool,
    p_mod_4: u8,
    /// `true` iff `−1` is a non-residue, i.e. `b ↦ −b` is the quadratic twist.
    sign_choice_selects_twist: bool,
    /// Which of the two roots of `b² ≡ a³/r` the standard published.
    published_root_is_smaller: bool,
    twist_order_is_prime: bool,
    twist_small_factors: Vec<(u64, u32)>,
    twist_residue_bits: u64,
    twist_residue_is_probable_prime: bool,
    /// log2 of the product of twist factors below the trial bound.
    twist_leak_bits_below_bound: f64,
    /// Twist order factored further with bounded Pollard rho
    /// (`pollard_iters` iterations per split); each entry is
    /// `(factor hex, exponent, probable prime?)`.  Exponent `0` is
    /// `deep_factor`'s marker for a residue the rho budget did not split.
    twist_deep_factors: Vec<(String, u32, bool)>,
    /// Bit length of the largest *proved-probable-prime* twist factor
    /// found; rho on the twist costs about half this many bits.
    twist_largest_prime_factor_bits: u64,
}

/// Slot S4: generator provenance.
#[derive(Serialize, Clone)]
struct GeneratorProvenance {
    generator_on_curve: bool,
    /// Smallest `x ≥ 0` with `x³ + ax + b` a square mod `p`.
    smallest_valid_x: String,
    generator_x_is_smallest_valid_x: bool,
    generator_y_is_smaller_root: bool,
    /// How many `x` values were scanned before the first valid one.
    x_values_scanned: u32,
}

/// Slot S6: smoothness of `n ± 1` (no known ECDLP consequence; recorded
/// so the slot has a measured row rather than an assumption).
#[derive(Serialize, Clone)]
struct OrderNeighbours {
    n_minus_1_small_factors: Vec<(u64, u32)>,
    n_minus_1_residue_bits: u64,
    n_minus_1_residue_is_probable_prime: bool,
    n_plus_1_small_factors: Vec<(u64, u32)>,
    n_plus_1_residue_bits: u64,
    n_plus_1_residue_is_probable_prime: bool,
}

fn order_neighbours(n: &BigUint, trial_bound: u64) -> OrderNeighbours {
    let (fm, rm) = trial_factor(&(n - BigUint::one()), trial_bound);
    let (fp, rp) = trial_factor(&(n + BigUint::one()), trial_bound);
    OrderNeighbours {
        n_minus_1_small_factors: fm,
        n_minus_1_residue_bits: rm.bits(),
        n_minus_1_residue_is_probable_prime: rm > BigUint::one()
            && crypto_lib::asymmetric::rsa::is_prime(&rm),
        n_plus_1_small_factors: fp,
        n_plus_1_residue_bits: rp.bits(),
        n_plus_1_residue_is_probable_prime: rp > BigUint::one()
            && crypto_lib::asymmetric::rsa::is_prime(&rp),
    }
}

/// Slot S7: is any published seed a SHA-1 image of a public string?
///
/// A fingerprint seed `SHA-1(secret ‖ counter)` is undetectable by
/// construction; this only checks the *public* nothing-up-my-sleeve
/// candidates a documented procedure would have used, plus relations
/// between the five seeds themselves.  Every match would be an
/// executed fact; every miss is one candidate family excluded.
#[derive(Serialize, Clone)]
struct SeedRelations {
    /// `(candidate description, matches the published seed)` pairs.
    checks: Vec<(String, bool)>,
    any_match: bool,
}

fn seed_relations(curve: &CurveParams, all: &[CurveParams]) -> SeedRelations {
    let seed = seed_for(curve.name);
    let mut checks: Vec<(String, bool)> = Vec::new();
    let mut push = |label: String, candidate: [u8; 20]| {
        checks.push((label, candidate == seed));
    };
    let names = [
        curve.name.to_string(),
        curve.name.to_lowercase(),
        curve.name.replace('-', ""),
        format!("secp{}r1", curve.p.bits()),
        format!("prime{}v1", curve.p.bits()),
        format!("NIST {}", curve.name),
        format!("{}", curve.p),
        format!("{:x}", curve.p),
        format!("{:X}", curve.p),
        format!("{}", curve.n),
        format!("{:x}", curve.n),
    ];
    for name in &names {
        push(format!("SHA-1({name:?})"), sha1(name.as_bytes()));
        push(
            format!("SHA-1({name:?} ‖ '\\n')"),
            sha1(format!("{name}\n").as_bytes()),
        );
    }
    for other in all {
        let other_seed = seed_for(other.name);
        if other.name != curve.name {
            push(format!("seed({})", other.name), other_seed);
        }
        push(format!("SHA-1(seed({}))", other.name), sha1(&other_seed));
        for delta in 1u8..=4 {
            let mut bumped = BigUint::from_bytes_be(&other_seed) + BigUint::from(delta);
            bumped %= BigUint::one() << 160u32;
            let mut bytes = [0u8; 20];
            let be = bumped.to_bytes_be();
            bytes[20 - be.len()..].copy_from_slice(&be);
            push(
                format!("SHA-1(seed({}) + {delta})", other.name),
                sha1(&bytes),
            );
            if other.name != curve.name {
                push(format!("seed({}) + {delta}", other.name), bytes);
            }
        }
    }
    for counter in 0u32..=16 {
        push(
            format!("SHA-1({:?} ‖ {counter} as ASCII)", curve.name),
            sha1(format!("{}{counter}", curve.name).as_bytes()),
        );
    }
    let any_match = checks.iter().any(|(_, m)| *m);
    SeedRelations { checks, any_match }
}

#[derive(Serialize, Clone)]
struct CurveCensus {
    name: String,
    p_bits: u64,
    trace_hex: String,
    trace_is_positive: bool,
    order_is_prime: bool,
    cofactor: u64,
    embedding_degree_if_at_most: u32,
    embedding_degree: Option<u32>,
    discriminant: DiscriminantProfile,
    root_sign: RootSignFreedom,
    generator: GeneratorProvenance,
    order_neighbours: OrderNeighbours,
    seed_relations: SeedRelations,
}

fn census(
    curve: &CurveParams,
    all: &[CurveParams],
    trial_bound: u64,
    pollard_iters: u64,
    max_k: u32,
) -> CurveCensus {
    let p = &curve.p;
    let report = curve_structural_report(curve, trial_bound);
    let seed = seed_for(curve.name);
    let derivation = verify_seed_derivation(curve, &seed);
    let (_, _, r) = x962_candidate_r(&seed, p.bits() as u32);

    // The two roots of b² ≡ a³ · r⁻¹.
    let a3 = (&curve.a * &curve.a % p) * &curve.a % p;
    let b_squared = (&a3 * mod_inverse(&r, p)) % p;
    let smaller_root = sqrt_mod(&b_squared, p).expect("r·b² = a³ must be solvable");
    let published_root_is_smaller = smaller_root == curve.b;
    let p_mod_4 = (p % 4u32).to_u8().unwrap();
    let minus_one = p - BigUint::one();
    let sign_choice_selects_twist = legendre(&minus_one, p) == -1;

    let twist_leak = {
        let mut acc = 0f64;
        for (q, e) in &report.twist_factors {
            acc += (*q as f64).log2() * (*e as f64);
        }
        acc
    };
    let twist_deep = deep_factor(&report.twist_order, trial_bound, pollard_iters);
    let twist_deep_factors: Vec<(String, u32, bool)> = twist_deep
        .iter()
        .map(|(q, e)| {
            (
                format!("0x{:x}", q),
                *e,
                crypto_lib::asymmetric::rsa::is_prime(q),
            )
        })
        .collect();
    let twist_largest_prime_factor_bits = twist_deep
        .iter()
        .filter(|(q, _)| crypto_lib::asymmetric::rsa::is_prime(q))
        .map(|(q, _)| q.bits())
        .max()
        .unwrap_or(0);
    let root_sign = RootSignFreedom {
        seed_derivation_holds: derivation.holds,
        p_mod_4,
        sign_choice_selects_twist,
        published_root_is_smaller,
        twist_order_is_prime: crypto_lib::asymmetric::rsa::is_prime(&report.twist_order),
        twist_small_factors: report.twist_factors.clone(),
        twist_residue_bits: report.twist_residue.bits(),
        twist_residue_is_probable_prime: report.twist_residue_is_prime,
        twist_leak_bits_below_bound: twist_leak,
        twist_deep_factors,
        twist_largest_prime_factor_bits,
    };

    // Generator provenance: is G the lexicographically first point?
    let rhs = |x: &BigUint| -> BigUint { ((x * x % p) * x % p + &curve.a * x % p + &curve.b) % p };
    let gy_on_curve = (&curve.gy * &curve.gy) % p == rhs(&curve.gx);
    let mut x = BigUint::zero();
    let mut scanned = 0u32;
    let smallest_valid_x = loop {
        scanned += 1;
        if sqrt_mod(&rhs(&x), p).is_some() {
            break x.clone();
        }
        x += BigUint::one();
    };
    let gy_root = sqrt_mod(&rhs(&curve.gx), p).expect("G is on the curve");
    let generator = GeneratorProvenance {
        generator_on_curve: gy_on_curve,
        smallest_valid_x: format!("0x{:x}", smallest_valid_x),
        generator_x_is_smallest_valid_x: smallest_valid_x == curve.gx,
        generator_y_is_smaller_root: gy_root == curve.gy,
        x_values_scanned: scanned,
    };

    let trace_is_positive = report.t > BigInt::zero();
    let embedding_degree = embedding_degree_bound(p, &curve.n, max_k);
    CurveCensus {
        name: curve.name.to_string(),
        p_bits: p.bits(),
        trace_hex: format!(
            "{}0x{:x}",
            if trace_is_positive { "" } else { "-" },
            report.t.magnitude()
        ),
        trace_is_positive,
        order_is_prime: crypto_lib::asymmetric::rsa::is_prime(&curve.n),
        cofactor: curve.h as u64,
        embedding_degree_if_at_most: max_k,
        embedding_degree,
        discriminant: discriminant_profile(&report.cm_disc_abs, trial_bound),
        root_sign,
        generator,
        order_neighbours: order_neighbours(&curve.n, trial_bound),
        seed_relations: seed_relations(curve, all),
    }
}

fn markdown(rows: &[CurveCensus]) -> String {
    let mut s = String::new();
    s.push_str("# NIST prime-curve trapdoor-slot census\n\n");
    s.push_str("Executed by `examples/nist_slot_census.rs`.  Every value is recomputed from `CurveParams` and the published seeds; nothing is pinned by hand.\n\n");
    s.push_str("## Group invariants (slots S2, S5, S6)\n\n");
    s.push_str("| curve | p bits | trace sign | n prime | k ≤ bound | \\|D\\| bits | square part of \\|D\\| from small factors (bits) | \\|D\\| residue (bits) | residue square? | residue prime? |\n|---|---:|---|---|---|---:|---:|---:|---|---|\n");
    for c in rows {
        s.push_str(&format!(
            "| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |\n",
            c.name,
            c.p_bits,
            if c.trace_is_positive { "+" } else { "−" },
            c.order_is_prime,
            c.embedding_degree
                .map(|k| k.to_string())
                .unwrap_or_else(|| format!("> {}", c.embedding_degree_if_at_most)),
            c.discriminant.abs_disc_bits,
            c.discriminant.square_part_bits_from_small_factors,
            c.discriminant.residue_bits,
            c.discriminant.residue_is_perfect_square,
            c.discriminant.residue_is_probable_prime,
        ));
    }
    s.push_str("\n## X9.62 root-sign freedom and twist (slots S3, S5)\n\n");
    s.push_str("| curve | seed → b holds | p mod 4 | sign choice = twist | published root smaller | twist order prime | twist leak bits (< trial bound) | twist residue bits | residue prime |\n|---|---|---:|---|---|---|---:|---:|---|\n");
    for c in rows {
        let r = &c.root_sign;
        s.push_str(&format!(
            "| {} | {} | {} | {} | {} | {} | {:.1} | {} | {} |\n",
            c.name,
            r.seed_derivation_holds,
            r.p_mod_4,
            r.sign_choice_selects_twist,
            r.published_root_is_smaller,
            r.twist_order_is_prime,
            r.twist_leak_bits_below_bound,
            r.twist_residue_bits,
            r.twist_residue_is_probable_prime
        ));
    }
    s.push_str("\n## Generator provenance (slot S4)\n\n");
    s.push_str("| curve | G on curve | smallest valid x | G.x is smallest valid x | G.y is smaller root | x scanned |\n|---|---|---|---|---|---:|\n");
    for c in rows {
        let g = &c.generator;
        s.push_str(&format!(
            "| {} | {} | {} | {} | {} | {} |\n",
            c.name,
            g.generator_on_curve,
            g.smallest_valid_x,
            g.generator_x_is_smallest_valid_x,
            g.generator_y_is_smaller_root,
            g.x_values_scanned
        ));
    }
    s.push_str("\n## Twist order factorisation (trial division + bounded Pollard rho)\n\n");
    s.push_str("| curve | factors (hex, exponent, probable prime) | largest proved prime factor (bits) |\n|---|---|---:|\n");
    for c in rows {
        let r = &c.root_sign;
        let fs: Vec<String> = r
            .twist_deep_factors
            .iter()
            .map(|(q, e, pr)| {
                if *e == 0 {
                    format!(
                        "{} (unsplit residue, {})",
                        q,
                        if *pr { "prime" } else { "composite" }
                    )
                } else {
                    format!("{}^{} ({})", q, e, if *pr { "prime" } else { "composite" })
                }
            })
            .collect();
        s.push_str(&format!(
            "| {} | {} | {} |\n",
            c.name,
            fs.join(" · "),
            r.twist_largest_prime_factor_bits
        ));
    }
    s.push_str("\n## |D| small factors\n\n");
    for c in rows {
        s.push_str(&format!(
            "- {}: {:?}\n",
            c.name, c.discriminant.small_factors
        ));
    }
    s.push_str("\n## n ± 1 smoothness (slot S6)\n\n");
    s.push_str("| curve | n−1 small factors | n−1 residue bits | residue prime | n+1 small factors | n+1 residue bits | residue prime |\n|---|---|---:|---|---|---:|---|\n");
    for c in rows {
        let o = &c.order_neighbours;
        s.push_str(&format!(
            "| {} | {:?} | {} | {} | {:?} | {} | {} |\n",
            c.name,
            o.n_minus_1_small_factors,
            o.n_minus_1_residue_bits,
            o.n_minus_1_residue_is_probable_prime,
            o.n_plus_1_small_factors,
            o.n_plus_1_residue_bits,
            o.n_plus_1_residue_is_probable_prime
        ));
    }
    s.push_str("\n## Seed relations (slot S7)\n\n");
    for c in rows {
        let r = &c.seed_relations;
        let hits: Vec<&String> = r
            .checks
            .iter()
            .filter(|(_, m)| *m)
            .map(|(l, _)| l)
            .collect();
        s.push_str(&format!(
            "- {}: {} candidates checked, {} matched{}\n",
            c.name,
            r.checks.len(),
            hits.len(),
            if hits.is_empty() {
                String::new()
            } else {
                format!(": {hits:?}")
            }
        ));
    }
    s
}

fn main() {
    let json = std::env::args().any(|a| a == "--json");
    let trial_bound: u64 = 1_000_000;
    let max_k: u32 = 1_000;
    let pollard_iters: u64 = 1 << 22;
    let curves = [
        CurveParams::p192(),
        CurveParams::p224(),
        CurveParams::p256(),
        CurveParams::p384(),
        CurveParams::p521(),
    ];
    let rows: Vec<CurveCensus> = curves
        .iter()
        .map(|c| census(c, &curves, trial_bound, pollard_iters, max_k))
        .collect();
    if json {
        println!("{}", serde_json::to_string_pretty(&rows).unwrap());
    } else {
        print!("{}", markdown(&rows));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn all() -> Vec<CurveParams> {
        vec![
            CurveParams::p192(),
            CurveParams::p224(),
            CurveParams::p256(),
            CurveParams::p384(),
            CurveParams::p521(),
        ]
    }

    #[test]
    fn x962_derivation_holds_for_every_published_seed() {
        for c in all() {
            assert!(
                verify_seed_derivation(&c, &seed_for(c.name)).holds,
                "{}",
                c.name
            );
        }
    }

    #[test]
    fn root_sign_selects_twist_exactly_when_p_is_3_mod_4() {
        let curves = all();
        for c in &curves {
            let row = census(c, &curves, 1_000, 1 << 10, 100);
            let expected = row.root_sign.p_mod_4 == 3;
            assert_eq!(
                row.root_sign.sign_choice_selects_twist, expected,
                "{}",
                c.name
            );
            assert_eq!(c.name == "P-224", row.root_sign.p_mod_4 == 1, "{}", c.name);
        }
    }

    #[test]
    fn published_root_choice_matches_round_zero() {
        let curves = all();
        let expected = [
            ("P-192", true),
            ("P-224", false),
            ("P-256", true),
            ("P-384", false),
            ("P-521", true),
        ];
        for (name, smaller) in expected {
            let c = curves.iter().find(|c| c.name == name).unwrap();
            let row = census(c, &curves, 1_000, 1 << 10, 100);
            assert_eq!(row.root_sign.published_root_is_smaller, smaller, "{name}");
        }
    }

    #[test]
    fn generator_is_never_the_first_valid_x() {
        let curves = all();
        for c in &curves {
            let row = census(c, &curves, 1_000, 1 << 10, 100);
            assert!(row.generator.generator_on_curve, "{}", c.name);
            assert!(!row.generator.generator_x_is_smallest_valid_x, "{}", c.name);
        }
    }

    #[test]
    fn tonelli_shanks_agrees_with_squaring() {
        let p = CurveParams::p224().p;
        let x = BigUint::from(123_456_789u64);
        let sq = (&x * &x) % &p;
        let r = sqrt_mod(&sq, &p).unwrap();
        assert_eq!((&r * &r) % &p, sq);
    }
}
