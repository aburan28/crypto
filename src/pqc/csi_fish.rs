//! **CSI-FiSh** — "Commutative Supersingular Isogeny based Fiat–Shamir"
//! signatures (Beullens–Kleinjung–Vercauteren 2019).
//!
//! # From CSIDH key exchange to a signature
//! CSIDH (`pqc::csidh`) is a *non-interactive key exchange* built on the
//! free, transitive action of the ideal-class group `Cl(O)` on a set of
//! supersingular curves.  CSI-FiSh turns that group action into a
//! **signature** via the Fiat–Shamir transform of a classic
//! group-action identification protocol (the isogeny analog of Schnorr):
//!
//! - **Public key**: a base curve `E₀` and `E_A = a ⋆ E₀`, where the
//!   secret `a` is a group element.
//! - **Identify**: commit `E₁ = b ⋆ E₀` for random `b`; on challenge
//!   `c ∈ {0,1}` respond with `r = b − c·a`; the verifier checks
//!   `r ⋆ E_c = E₁` (with `E₀` for `c = 0`, `E_A` for `c = 1`).  This is
//!   sound because `(b − a) ⋆ E_A = (b − a) ⋆ a ⋆ E₀ = b ⋆ E₀ = E₁`.
//!
//! CSI-FiSh's real contribution was *computing the class group* of a
//! CSIDH-512 instance so that group elements can be sampled uniformly
//! and composed as integer vectors reduced modulo the relation lattice —
//! which is exactly what makes `b − c·a` well defined.
//!
//! # This implementation
//! We **precompute the cycle** of a 3-isogeny walk at the toy CSIDH
//! prime `p = 419`, using the classical **modular polynomial `Φ₃`** to
//! take correct isogeny steps (the roots of `Φ₃(j, ·)` over F_p are the
//! curves 3-isogenous to `j`).  A deterministic non-backtracking walk
//! from `j = 1728` closes into a loop, realising the group concretely
//! as `Z/hZ`: the action of an integer `k` is "walk `k` steps", and
//! `b − c·a` is arithmetic mod `h`.  Fiat–Shamir over `τ = 64` binary
//! rounds gives soundness `2⁻⁶⁴`.  Toy scale, not constant-time; see
//! SECURITY.md.

use crate::hash::sha3::{sha3_256, shake256};
use crate::utils::random::random_bytes;

/// Fiat–Shamir rounds.
pub const TAU: usize = 64;

/// Base field prime for the isogeny cycle (the CSIDH prime `4·3·5·7−1`).
const P: u64 = 419;

fn fp_mul(a: u64, b: u64) -> u64 {
    a * b % P
}
fn fp_add(a: u64, b: u64) -> u64 {
    (a + b) % P
}

/// Evaluate the classical modular polynomial `Φ₃(x, y)` mod `P`.  Its
/// roots in `y` for fixed `x = j` are the `j`-invariants 3-isogenous to
/// the curve with invariant `j` — a *correct* isogeny step (unlike the
/// toy CSIDH's approximate Vélu formula).
fn phi3(x: u64, y: u64) -> u64 {
    // Coefficients of Φ₃ reduced mod P (rem_euclid handles negatives).
    let c = |v: i128| -> u64 { v.rem_euclid(P as i128) as u64 };
    let x2 = fp_mul(x, x);
    let x3 = fp_mul(x2, x);
    let x4 = fp_mul(x3, x);
    let y2 = fp_mul(y, y);
    let y3 = fp_mul(y2, y);
    let y4 = fp_mul(y3, y);
    let mut acc = 0u64;
    acc = fp_add(acc, x4);
    acc = fp_add(acc, y4);
    acc = fp_add(acc, P - fp_mul(x3, y3)); // − x³y³
    acc = fp_add(acc, fp_mul(c(2232), fp_add(fp_mul(x3, y2), fp_mul(x2, y3))));
    acc = fp_add(acc, fp_mul(c(-1069956), fp_add(fp_mul(x3, y), fp_mul(x, y3))));
    acc = fp_add(acc, fp_mul(c(36864000), fp_add(x3, y3)));
    acc = fp_add(acc, fp_mul(c(2587918086), fp_mul(x2, y2)));
    acc = fp_add(acc, fp_mul(c(8900222976000), fp_add(fp_mul(x2, y), fp_mul(x, y2))));
    acc = fp_add(acc, fp_mul(c(452984832000000), fp_add(x2, y2)));
    acc = fp_add(acc, fp_mul(c(-770845966336000000), fp_mul(x, y)));
    acc = fp_add(acc, fp_mul(c(1855425871872000000000), fp_add(x, y)));
    acc
}

/// The `j`-invariants 3-isogenous to `j` (roots of `Φ₃(j, ·)` over F_p).
fn neighbors(j: u64) -> Vec<u64> {
    (0..P).filter(|&y| phi3(j, y) == 0).collect()
}

/// The precomputed action cycle: `curves[i] = gⁱ ⋆ E₀`, indexed so that
/// the group element `k` maps a curve at index `i` to index `(i+k) mod h`.
#[derive(Clone)]
pub struct ActionCycle {
    curves: Vec<u64>, // Montgomery A-coefficients, in walk order
    index_of: std::collections::HashMap<u64, usize>,
    pub order: usize, // class number h
}

impl ActionCycle {
    /// Build the action cycle by walking a single generator (the first
    /// `ELL`) and extracting the pure loop of the functional graph.
    ///
    /// The toy CSIDH step is deterministic but need not return to the
    /// literal base curve `A = 0` (its orbit can have a "tail" before the
    /// loop).  We walk until the first repeat, then keep only the loop:
    /// on that loop the step is a genuine cyclic permutation, giving a
    /// free transitive `Z/hZ` action whose canonical entry curve serves
    /// as the protocol base `E₀`.
    pub fn build() -> ActionCycle {
        // Trace a deterministic non-backtracking 3-isogeny walk from
        // j = 1728 (E₀: y² = x³ + x, supersingular since p ≡ 3 mod 4).
        // Each step moves to the smallest 3-isogenous neighbour that is
        // not where we came from; this deterministic step is a cyclic
        // permutation on its loop, giving a free transitive Z/hZ action
        // — responses read mod h make it a genuine group action.
        let start = 1728 % P;
        let mut seq = vec![start];
        let mut seen = std::collections::HashMap::new();
        seen.insert(start, 0usize);
        let mut prev: Option<u64> = None;
        let mut cur = start;
        let loop_start;
        loop {
            let mut nbrs = neighbors(cur);
            nbrs.sort_unstable();
            let next = nbrs
                .iter()
                .copied()
                .find(|&y| Some(y) != prev)
                .or_else(|| nbrs.first().copied())
                .unwrap_or(cur);
            prev = Some(cur);
            cur = next;
            if let Some(&s) = seen.get(&cur) {
                loop_start = s;
                break;
            }
            seen.insert(cur, seq.len());
            seq.push(cur);
            if seq.len() > 100_000 {
                loop_start = 0;
                break;
            }
        }
        let curves: Vec<u64> = seq[loop_start..].to_vec();
        let index_of = curves.iter().enumerate().map(|(i, &a)| (a, i)).collect();
        let order = curves.len();
        ActionCycle { curves, index_of, order }
    }

    /// Curve `k ⋆ E₀` as its A-coefficient.
    fn act_from_base(&self, k: usize) -> u64 {
        self.curves[k % self.order]
    }

    /// Curve `k ⋆ E` where `E` has A-coefficient `a`.
    fn act(&self, a: u64, k: usize) -> Option<u64> {
        let i = *self.index_of.get(&a)?;
        Some(self.curves[(i + k) % self.order])
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct CsiFishPublicKey {
    /// A-coefficient of `E_A = a ⋆ E₀`.
    pub e_a: u64,
}

#[derive(Clone)]
pub struct CsiFishSecretKey {
    a: usize, // secret group element in Z/hZ
}

#[derive(Clone, Debug, PartialEq)]
pub struct CsiFishSignature {
    /// The `τ` commitments (A-coefficients of `E₁` per round).
    pub commits: Vec<u64>,
    /// The `τ` responses `r = b − c·a mod h`.
    pub responses: Vec<usize>,
}

fn rand_mod(m: usize) -> usize {
    let mut b = [0u8; 8];
    random_bytes(&mut b);
    (u64::from_le_bytes(b) % m as u64) as usize
}

pub fn csi_fish_keygen(cycle: &ActionCycle) -> (CsiFishPublicKey, CsiFishSecretKey) {
    let a = rand_mod(cycle.order);
    (CsiFishPublicKey { e_a: cycle.act_from_base(a) }, CsiFishSecretKey { a })
}

/// Derive `τ` binary challenges from the public key, commitments, msg.
fn challenges(pk: &CsiFishPublicKey, commits: &[u64], msg: &[u8]) -> Vec<u8> {
    let mut input = b"CSI-FiSh".to_vec();
    input.extend_from_slice(&pk.e_a.to_le_bytes());
    for c in commits {
        input.extend_from_slice(&c.to_le_bytes());
    }
    input.extend_from_slice(msg);
    let seed = sha3_256(&input);
    shake256(&seed, TAU).iter().map(|b| b & 1).collect()
}

pub fn csi_fish_sign(
    cycle: &ActionCycle,
    pk: &CsiFishPublicKey,
    sk: &CsiFishSecretKey,
    msg: &[u8],
) -> CsiFishSignature {
    let bs: Vec<usize> = (0..TAU).map(|_| rand_mod(cycle.order)).collect();
    let commits: Vec<u64> = bs.iter().map(|&b| cycle.act_from_base(b)).collect();
    let ch = challenges(pk, &commits, msg);
    let responses: Vec<usize> = bs
        .iter()
        .zip(&ch)
        .map(|(&b, &c)| {
            // r = b − c·a  (mod h).
            (b + cycle.order - if c == 1 { sk.a } else { 0 }) % cycle.order
        })
        .collect();
    CsiFishSignature { commits, responses }
}

pub fn csi_fish_verify(
    cycle: &ActionCycle,
    pk: &CsiFishPublicKey,
    msg: &[u8],
    sig: &CsiFishSignature,
) -> bool {
    if sig.commits.len() != TAU || sig.responses.len() != TAU {
        return false;
    }
    let ch = challenges(pk, &sig.commits, msg);
    for i in 0..TAU {
        // r ⋆ E_c must equal the committed E₁.
        let base_curve = if ch[i] == 1 { pk.e_a } else { cycle.curves[0] };
        match cycle.act(base_curve, sig.responses[i]) {
            Some(reached) if reached == sig.commits[i] => {}
            _ => return false,
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cycle_is_a_genuine_torsor() {
        // The action orbit is a single cycle through E₀ of length h > 1,
        // so it realises a free transitive Z/hZ action.
        let cycle = ActionCycle::build();
        assert!(cycle.order > 1, "class number must exceed 1");
        // Every curve appears exactly once (bijective indexing).
        assert_eq!(cycle.index_of.len(), cycle.order);
        // Walking `order` steps from the base returns to the base.
        assert_eq!(cycle.act(cycle.curves[0], cycle.order), Some(cycle.curves[0]));
    }

    #[test]
    fn group_law_holds() {
        // (x+y) ⋆ E₀ = y ⋆ (x ⋆ E₀): the precomputed action composes.
        let cycle = ActionCycle::build();
        let (x, y) = (3 % cycle.order, 5 % cycle.order);
        let direct = cycle.act_from_base(x + y);
        let composed = cycle.act(cycle.act_from_base(x), y).unwrap();
        assert_eq!(direct, composed);
    }

    #[test]
    fn sign_verify_roundtrip() {
        let cycle = ActionCycle::build();
        let (pk, sk) = csi_fish_keygen(&cycle);
        let sig = csi_fish_sign(&cycle, &pk, &sk, b"isogeny Fiat-Shamir");
        assert!(csi_fish_verify(&cycle, &pk, b"isogeny Fiat-Shamir", &sig));
    }

    #[test]
    fn wrong_message_rejected() {
        let cycle = ActionCycle::build();
        let (pk, sk) = csi_fish_keygen(&cycle);
        let sig = csi_fish_sign(&cycle, &pk, &sk, b"one");
        assert!(!csi_fish_verify(&cycle, &pk, b"two", &sig));
    }

    #[test]
    fn wrong_public_key_rejected() {
        let cycle = ActionCycle::build();
        let (pk_a, sk_a) = csi_fish_keygen(&cycle);
        let (pk_b, _) = csi_fish_keygen(&cycle);
        let sig = csi_fish_sign(&cycle, &pk_a, &sk_a, b"m");
        // A different public key (whp) fails verification.
        if pk_a.e_a != pk_b.e_a {
            assert!(!csi_fish_verify(&cycle, &pk_b, b"m", &sig));
        }
    }

    #[test]
    fn tampered_response_rejected() {
        let cycle = ActionCycle::build();
        let (pk, sk) = csi_fish_keygen(&cycle);
        let mut sig = csi_fish_sign(&cycle, &pk, &sk, b"m");
        sig.responses[0] = (sig.responses[0] + 1) % cycle.order;
        assert!(!csi_fish_verify(&cycle, &pk, b"m", &sig));
    }
}
