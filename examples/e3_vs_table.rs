//! Lead 3 of the prime-degree survey (`research/pdp_shape_replay_20261015`,
//! `EXPERIMENTS.md`): the exact `m = 3` decomposition by enumerating `e₃`
//! (`koblitz_sym_s4::decompose_by_e3`), counted in `F_{2^n}`
//! multiplications per target, against the root-finding pair-table
//! oracle's `5|F| − 3` multiplications per target with `|F| = 2^{ℓ−1}`
//! (`research/pdp_root_finder_20261012` §0, the same unit).  Koblitz
//! curve `b = 1`, `V = ⟨1, z, …, z^{ℓ−1}⟩` (geometric: `θ = 1`, `g = z`),
//! uniform random target abscissae.  One JSON line per cell on stdout.
//!
//! ```sh
//! cargo run --release --example e3_vs_table -- 16
//! ```

use std::time::Instant;

use crypto_lib::binary_ecc::F2mElement;
use crypto_lib::cryptanalysis::koblitz_index_calculus::find_irreducible_sparse;
use crypto_lib::cryptanalysis::koblitz_sym_s4::{decompose_by_e3, SymSpaces, WordField};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::json;

fn main() {
    let targets: u32 = std::env::args()
        .nth(1)
        .map_or(16, |s| s.parse().expect("targets"));
    for (n, l) in [
        (13u32, 4u32),
        (13, 5),
        (17, 4),
        (17, 5),
        (19, 5),
        (19, 6),
        (23, 5),
        (23, 6),
    ] {
        let irr = find_irreducible_sparse(n).expect("an irreducible");
        let f = WordField::new(n, &irr);
        let one = F2mElement::one(n);
        let z = F2mElement::from_bit_positions(&[1], n);
        let spaces = SymSpaces::geometric(n, l, &one, &z, &irr);
        let v: Vec<u64> = (0..1u64 << l).collect();
        let mut rng = StdRng::seed_from_u64(20261015 ^ (u64::from(n) << 8) ^ u64::from(l));
        let (mut muls, mut xors, mut visited, mut roots, mut split, mut found) =
            (0u64, 0u64, 0u64, 0u64, 0u64, 0u32);
        let t = Instant::now();
        for _ in 0..targets {
            let x_r = rng.gen::<u64>() & ((1u64 << n) - 1);
            let (hit, st) = decompose_by_e3(&spaces, &v, x_r, 1, &f, |_| true);
            muls += st.field_muls;
            xors += st.word_xors;
            visited += st.e3_visited;
            roots += st.roots;
            split += st.split;
            found += u32::from(hit);
        }
        let wall = t.elapsed().as_secs_f64() / f64::from(targets);
        let table = 5 * (1u64 << (l - 1)) - 3;
        let per = muls as f64 / f64::from(targets);
        println!(
            "{}",
            json!({
                "n": n, "l": l, "targets": targets, "decomposed": found,
                "dim_w": spaces.w.iter().map(|w| w.len()).collect::<Vec<_>>(),
                "e3_visited_per_target": visited as f64 / f64::from(targets),
                "e3_field_muls_per_target": per,
                "e3_word_xors_per_target": xors as f64 / f64::from(targets),
                "roots": roots, "split": split,
                "table_muls_per_target": table,
                "e3_over_table": per / table as f64,
                "e3_wall_s_per_target": wall,
            })
        );
    }
}
