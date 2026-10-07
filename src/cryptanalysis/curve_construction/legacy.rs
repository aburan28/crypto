//! Native replay of the numbers `research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md`
//! took from its legacy Python script: curve constants and the rho reference,
//! the trace-zero polynomial with `#A(F_2) = r`, the GHS magic-number
//! trichotomy (through the library's `ec_trapdoor` magic number), the trace
//! annihilating `⟨G⟩`, and the zeta-function window.

use super::weil::{
    curve_order, eval, frobenius_trace, is_probable_prime, log2_big, sqrt_mod, trace_zero_poly,
};
use super::window::{price, rational_model_points, Cell};
use crate::binary_ecc::{F2mElement, IrreduciblePoly};
use crate::cryptanalysis::ec_trapdoor::{ghs_genus_with_type, magic_number_full, FieldTower};
use num_bigint::{BigInt, BigUint};
use num_traits::{One, Zero};
use rand::{rngs::StdRng, Rng, SeedableRng};
use serde::Serialize;
use std::collections::BTreeMap;

pub const N: u32 = 131;
pub const AUTOMORPHISMS: u64 = 2 * N as u64;

#[derive(Clone, Debug, Serialize)]
pub struct ChallengeTarget {
    pub curve: String,
    pub order: String,
    pub cofactor: u32,
    pub r: String,
    pub r_is_prime: bool,
    pub log2_r: f64,
    pub rho_automorphism_order: u64,
    pub log2_rho_plain: f64,
    pub log2_rho_reference: f64,
    pub s_rho_reference: f64,
}

pub fn challenge_target() -> (ChallengeTarget, BigUint) {
    let order = curve_order(-1, N).to_biguint().expect("positive");
    let r: BigUint = &order / 4u32;
    assert!((&order % 4u32).is_zero());
    let log2r = log2_big(&r);
    let log2_pi = std::f64::consts::PI.log2();
    let plain = 0.5 * (log2_pi + log2r - 1.0);
    let reference = 0.5 * (log2_pi + log2r - (2.0 * AUTOMORPHISMS as f64).log2());
    (
        ChallengeTarget {
            curve: "y^2 + x y = x^3 + 1 over F_2, used over F_2^131 (ECC2K-130)".into(),
            order: order.to_string(),
            cofactor: 4,
            r: r.to_string(),
            r_is_prime: is_probable_prime(&r),
            log2_r: log2r,
            rho_automorphism_order: AUTOMORPHISMS,
            log2_rho_plain: plain,
            log2_rho_reference: reference,
            s_rho_reference: 2f64.powf(reference - log2r / 2.0),
        },
        r,
    )
}

#[derive(Clone, Debug, Serialize)]
pub struct TraceZeroFacts {
    pub s_131: String,
    pub degree: usize,
    pub points_on_a_equal_r: bool,
    pub functional_equation_holds: bool,
    pub conductor_7_divides_131: bool,
}

#[derive(Clone, Debug, Serialize)]
pub struct GhsCensus {
    pub field_modulus: String,
    pub ecc2k130_magic_number: u32,
    pub ecc2k130_descent_genus: String,
    pub seed: u64,
    pub samples: usize,
    pub census: BTreeMap<String, u64>,
    pub genus_by_magic_number: BTreeMap<u32, String>,
    pub trichotomy_holds: bool,
}

#[derive(Clone, Debug, Serialize)]
pub struct TraceKill {
    pub lambda: String,
    pub lambda_satisfies_char_poly: bool,
    pub lambda_pow_131_is_one: bool,
    pub sum_lambda_i_is_zero_mod_r: bool,
}

#[derive(Clone, Debug, Serialize)]
pub struct WindowReplay {
    pub exact_genus_130_jac_a: Cell,
    pub cells: Vec<Cell>,
    pub crossover_between: (usize, usize),
}

#[derive(Clone, Debug, Serialize)]
pub struct LegacyReplay {
    pub target: ChallengeTarget,
    pub trace_zero: TraceZeroFacts,
    pub ghs: GhsCensus,
    pub trace_kill: TraceKill,
    pub window: WindowReplay,
}

fn f2_trace(b: &F2mElement, irr: &IrreduciblePoly) -> F2mElement {
    let mut acc = b.clone();
    let mut cur = b.clone();
    for _ in 1..N {
        cur = cur.square(irr);
        acc = acc.add(&cur);
    }
    acc
}

pub fn ghs_census(samples: usize, seed: u64) -> GhsCensus {
    let irr = IrreduciblePoly::deg_131();
    let tower = FieldTower::new(N, N, 1, irr.clone());
    let zero = F2mElement::zero(N);
    let one = F2mElement::one(N);
    let m_challenge = magic_number_full(&tower, &zero, &one);
    let (g_challenge, _) = ghs_genus_with_type(&tower, &one);
    let mut rng = StdRng::seed_from_u64(seed);
    let mut census = BTreeMap::new();
    let mut genus_by_magic = BTreeMap::new();
    let mut ok = true;
    for _ in 0..samples {
        let v = loop {
            let lo: u128 = rng.gen();
            let hi: u8 = rng.gen::<u8>() & 7;
            let v: BigUint = (BigUint::from(hi) << 128u32) | BigUint::from(lo);
            if !v.is_zero() {
                break v;
            }
        };
        let b = F2mElement::from_biguint(&v, N);
        let m = magic_number_full(&tower, &zero, &b);
        let tr = f2_trace(&b, &irr);
        assert!(tr == zero || tr == one, "absolute trace lies in F_2");
        let tr_bit = u8::from(tr == one);
        ok &= matches!(m, 1 | 130 | 131);
        ok &= (m == 1) == (b == one);
        ok &= (m == 130) == (tr_bit == 0 && b != one);
        *census.entry(format!("m={m},trace={tr_bit}")).or_insert(0) += 1;
        genus_by_magic.entry(m).or_insert_with(|| {
            let (g, type_i) = ghs_genus_with_type(&tower, &b);
            format!("{g} (type {})", if type_i { "I" } else { "II" })
        });
    }
    GhsCensus {
        field_modulus: "z^131 + z^8 + z^3 + z^2 + 1 (IrreduciblePoly::deg_131)".into(),
        ecc2k130_magic_number: m_challenge,
        ecc2k130_descent_genus: g_challenge.to_string(),
        seed,
        samples,
        census,
        genus_by_magic_number: genus_by_magic,
        trichotomy_holds: ok,
    }
}

pub fn trace_kill(r: &BigUint) -> TraceKill {
    let root = sqrt_mod(&(r - 7u32), r).expect("−7 is a square mod r: E_0 has CM by Q(√−7)");
    let inv2 = (r + 1u32) >> 1;
    let rm1: BigUint = r - 1u32;
    let cands: [BigUint; 2] = [(&rm1 + &root) * &inv2 % r, (&rm1 + (r - &root)) * &inv2 % r];
    let lam = cands
        .iter()
        .find(|l| l.modpow(&BigUint::from(N), r).is_one())
        .expect("one eigenvalue has order dividing 131")
        .clone();
    let char_ok = ((&lam * &lam + &lam + 2u32) % r).is_zero();
    let mut sum = BigUint::zero();
    let mut pw = BigUint::one();
    for _ in 0..N {
        sum = (sum + &pw) % r;
        pw = pw * &lam % r;
    }
    TraceKill {
        lambda: lam.to_string(),
        lambda_satisfies_char_poly: char_ok,
        lambda_pow_131_is_one: lam.modpow(&BigUint::from(N), r).is_one(),
        sum_lambda_i_is_zero_mod_r: sum.is_zero(),
    }
}

/// `#C(F_{2^k})`, `k ≤ 130`, for a genus-130 curve with `Jac(C) ~ A`:
/// `tr_k(A) = −s_k` for `131 ∤ k`.
pub fn a_points() -> Vec<BigInt> {
    let mut v = vec![BigInt::zero()];
    v.extend((1..=130u32).map(|k| (BigInt::one() << k) + 1 + frobenius_trace(-1, k)));
    v
}

pub const WINDOW_GENERA: [usize; 13] = [
    130, 140, 150, 175, 200, 225, 250, 275, 290, 300, 325, 350, 400,
];

pub fn run(samples: usize, seed: u64) -> LegacyReplay {
    let (target, r) = challenge_target();
    let pa = trace_zero_poly(-1, N);
    let pa_at_1 = eval(&pa, &BigInt::one());
    let fe = (0..=130usize).all(|i| pa[i] == (&pa[260 - i] << (130 - i)));
    let trace_zero = TraceZeroFacts {
        s_131: frobenius_trace(-1, N).to_string(),
        degree: pa.len() - 1,
        points_on_a_equal_r: pa_at_1 == BigInt::from(r.clone()),
        functional_equation_holds: fe,
        conductor_7_divides_131: N.is_multiple_of(7),
    };
    let rho = target.log2_rho_reference;
    let exact = price(
        &a_points(),
        130,
        44,
        "exact: place counts from the zeta function of A",
        rho,
    );
    let cells: Vec<Cell> = WINDOW_GENERA
        .iter()
        .map(|&g| price(&rational_model_points(g), g, 44, "random-polynomial", rho))
        .collect();
    let last_win = cells
        .iter()
        .filter(|c| c.beats_rho)
        .map(|c| c.genus)
        .max()
        .expect("genus 130 wins");
    let first_loss = cells
        .iter()
        .filter(|c| !c.beats_rho)
        .map(|c| c.genus)
        .min()
        .expect("genus 400 loses");
    LegacyReplay {
        ghs: ghs_census(samples, seed),
        trace_kill: trace_kill(&r),
        target,
        trace_zero,
        window: WindowReplay {
            exact_genus_130_jac_a: exact,
            cells,
            crossover_between: (last_win, first_loss),
        },
    }
}
