//! Steps 2–3 of the fixed-degree algorithm: a connecting ideal between two maximal orders and
//! its normalised norm form.
//!
//! * `connecting_ideal` (Kirschmer–Voight, "Algorithmic enumeration of ideal classes for
//!   quaternion orders", SIAM J. Comput. 2010, §1–3): for maximal orders O_1, O_2 of B_{p,inf}
//!   the lattice O_1 O_2 is a fractional ideal with left order O_1 and right order O_2 and
//!   reduced norm 1/N, where N is the index [O_1 : O_1 ∩ O_2]; I = N O_1 O_2 is then an integral
//!   connecting ideal of norm N. (N is read off the reduced norm of O_1 O_2 rather than computed
//!   as an index; `connects` checks both orders.)
//! * `connecting_via_o0`: the same connection through E_0, conj(I_1) I_2 for left O_0-ideals I_1,
//!   I_2 (used as an independent cross-check: both have left order O_R(I_1) and right order
//!   O_R(I_2)).
//! * `norm_form`: an LLL-reduced basis sigma_1..4 of an integral left O_1-ideal I and the integral
//!   Gram matrix G of 2 Nrd(x) / N(I) on it (the slides' g_ij = tr(sigma_i conj(sigma_j)) / N(I)),
//!   so Q(x) = sum_ij G_ij x_i x_j = 2 Nrd(sum x_i sigma_i) / N(I).
use crate::int::Int;
use crate::quat::{Alg, Lattice, Quat};

/// I = N O_1 O_2 with Nrd(O_1 O_2) = 1/N (Kirschmer–Voight).
pub fn connecting_ideal(alg: &Alg, o1: &Lattice, o2: &Lattice) -> Lattice {
    let l = o1.mul(alg, o2);
    let (a, b) = l.norm(o1);
    // N(l) = a / b; scaling by t multiplies the norm by t^2: N(t l) = t^2 a / b = N  when t = b / a
    let i = l.scale(&b, &a);
    let (n, d) = i.norm(o1);
    assert!(d == Int::one(), "connecting ideal not integral: norm {n}/{d}");
    i
}

/// conj(I_1) I_2 for left O_0-ideals I_1, I_2 (left order O_R(I_1), right order O_R(I_2)).
pub fn connecting_via_o0(alg: &Alg, i1: &Lattice, i2: &Lattice) -> Lattice {
    i1.conj().mul(alg, i2)
}

/// Is `i` an ideal with exactly these left and right orders?
pub fn connects(alg: &Alg, i: &Lattice, o1: &Lattice, o2: &Lattice) -> bool {
    &i.left_order(alg, o1) == o1 && &i.right_order(alg, o1) == o2
}

#[derive(Clone, Debug)]
pub struct NormForm {
    /// LLL-reduced basis of I (quaternions)
    pub basis: Vec<Quat>,
    /// Q(x) = sum G_ij x_i x_j = 2 Nrd(sum x_i b_i) / N(I), integral and even on the diagonal
    pub gram: [[Int; 4]; 4],
    /// N(I)
    pub norm: Int,
}

/// The normalised norm form of an integral left O-ideal I (Step 3).
pub fn norm_form(alg: &Alg, o_left: &Lattice, i: &Lattice) -> NormForm {
    let (n, d) = i.norm(o_left);
    assert!(d == Int::one(), "norm_form: integral ideal expected");
    let rb = i.reduced_basis(alg);
    let den2 = &i.den * &i.den;
    let scale = &den2 * &n;
    let two = Int::from(2i64);
    let gram: [[Int; 4]; 4] = std::array::from_fn(|r| {
        std::array::from_fn(|s| {
            let num = &two * &alg.bil(&rb[r], &rb[s]);
            let (q, rem) = num.divrem_trunc(&scale);
            assert!(rem.is_zero(), "norm form not integral");
            q
        })
    });
    NormForm {
        basis: rb.iter().map(|v| Quat::new(v.clone(), i.den.clone())).collect(),
        gram,
        norm: n,
    }
}

/// sum x_i b_i.
pub fn combine(basis: &[Quat], x: &[Int; 4]) -> Quat {
    let mut a = Quat::zero();
    for (b, c) in basis.iter().zip(x) {
        if !c.is_zero() {
            a = a.add(&b.scale(c, &Int::one()));
        }
    }
    a
}

/// The normalised form of a maximal order O itself (Nrd on O: Q(x) = 2 Nrd(x)), for finding
/// endomorphisms of a prescribed degree.
pub fn order_norm_form(alg: &Alg, o: &Lattice) -> NormForm {
    norm_form(alg, o, o)
}
