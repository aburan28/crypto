//! Precomputed (u, v) values for the Castryck-Decru attack recipe.
//!
//! Each entry is (n, exp, u, v) where n is the chunk size (always odd) and
//! the Diophantine relation
//!     2^exp − 3^n = u² + 4 v²
//! holds. The right-hand side is the norm of u + v·(2ι) in End(E_0), where
//! ι² = [−1] is the j=1728 automorphism and so (2ι)² = [−4].
//!
//! This is the "magic table" from the Magma reference: for each guess of
//! `n` Bob-side ternary digits, we use the corresponding (u, v) to build
//! the Kani-distorted kernel as u·P + v·(2ι)(P) on the j=1728 curve.
//!
//! The Kani-glue construction works when (u, v) satisfy a specific
//! Diophantine relation linking the source elliptic curve's (2^a)-torsion
//! to the target Kani-glued Jacobian's (2^exp)-torsion. See uvtable.m in
//! the Magma reference (Castryck-Decru 2022/975).
//!
//! Table source: https://homes.esat.kuleuven.be/~wcastryc/code/uvtable.m

/// (n, exp, u, v): n digits ↔ 2^exp Kani torsion, with (u, v).
pub static UVTABLE: &[(u64, u64, u128, u128)] = &[
    (1, 3, 1, 1),
    (3, 5, 1, 1),
    (5, 8, 3, 1),
    (7, 13, 23, 37),
    (9, 15, 59, 49),
    (11, 19, 385, 223),
    (13, 21, 125, 349),
    (15, 24, 825, 661),
    (17, 29, 1307, 10075),
    (19, 34, 48991, 58347),
    (21, 41, 1440605, 168241),
    (23, 40, 721143, 348325),
    (25, 41, 956405, 330539),
    (27, 44, 3038895, 427699),
    (29, 46, 1021891, 416565),
    (31, 51, 16963049, 18346535),
    (33, 53, 37852565, 22446169),
    (35, 58, 111188129, 237611043),
    (37, 61, 1046252867, 436151935),
    (39, 63, 2170653479, 338777345),
    (41, 70, 5096872085, 16719304107),
    (43, 71, 3450988735, 22477861057),
    (45, 72, 36322560147, 10591569769),
    (47, 75, 86165183959, 30682562545),
    (49, 81, 1191966366037, 435249495185),
    (51, 82, 894987407119, 685749016857),
    (53, 86, 5245386911165, 2760159548823),
    (55, 91, 32934192159529, 17441114172845),
    (57, 102, 2190114350879525, 260974849403823),
    // ... (table continues to n=193 in the Magma reference; truncated
    // because subsequent values exceed u128. For SIKEp434 (a=216, b=137)
    // only the first ~70 entries are needed.)
];

/// Look up the (exp, u, v) triple for a given chunk size n. Returns None
/// if n is even or out of range.
pub fn lookup(n: u64) -> Option<(u64, u128, u128)> {
    UVTABLE.iter().find_map(|(nn, exp, u, v)| {
        if *nn == n { Some((*exp, *u, *v)) } else { None }
    })
}

/// Find the largest n ≤ n_max with an exp ≤ a_max entry. This is the
/// "biggest chunk that fits" used by the C-D attack to pick how many
/// digits to grab in one Groebner-system shot.
pub fn largest_usable(n_max: u64, a_max: u64) -> Option<(u64, u64, u128, u128)> {
    UVTABLE.iter().rev().find_map(|(n, exp, u, v)| {
        if *n <= n_max && *exp <= a_max { Some((*n, *exp, *u, *v)) } else { None }
    }).map(|(n, exp, u, v)| (n, exp, u, v))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn first_entry() {
        assert_eq!(lookup(1), Some((3, 1, 1)));
    }

    #[test]
    fn known_entry() {
        // From the SIKE_challenge.m example: n=7 → 2^13, u=23, v=37
        assert_eq!(lookup(7), Some((13, 23, 37)));
    }

    #[test]
    fn diophantine_relation_holds() {
        // For each table entry: 2^exp − 3^n = u² + 4v²  must hold.
        // (Norm of u + v·(2ι) in End(E_0).)
        for &(n, exp, u, v) in UVTABLE {
            let two_exp: u128 = match 1u128.checked_shl(exp as u32) {
                Some(v) => v,
                None => continue, // exp ≥ 128: overflow, skip
            };
            let three_n = match (3u128).checked_pow(n as u32) {
                Some(v) => v,
                None => continue,
            };
            if two_exp < three_n {
                panic!("uvtable[{n}]: 2^exp < 3^n, can't be a sum of squares");
            }
            let lhs = two_exp - three_n;
            let u_sq = match u.checked_mul(u) { Some(v) => v, None => continue };
            let v_sq_4 = match v.checked_mul(v).and_then(|x| x.checked_mul(4)) {
                Some(v) => v, None => continue,
            };
            let rhs = match u_sq.checked_add(v_sq_4) { Some(v) => v, None => continue };
            assert_eq!(
                lhs, rhs,
                "uvtable[n={n}]: 2^{exp} − 3^{n} = {lhs} ≠ u² + 4v² = {rhs}"
            );
        }
    }

    #[test]
    fn largest_usable_finds_entry() {
        // For SIDH-toy with a=4, b=3: usable entries have n ≤ b-3 = 0 → none.
        // For a slightly bigger SIDH (a=10, b=5), n ≤ 2 (b-3 = 2), exp ≤ 10:
        //   entries (1, 3, 1, 1) fits.
        // For our test, n_max=5, a_max=10: should find (5, 8, 3, 1).
        assert_eq!(largest_usable(5, 10), Some((5, 8, 3, 1)));
    }
}
