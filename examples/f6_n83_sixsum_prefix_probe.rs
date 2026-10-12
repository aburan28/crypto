//! Frozen one-prefix gate: research/f6_n83_sixsum_prefix_20261005.
use std::collections::HashSet;
use std::io::Write;
use std::time::Instant;

use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_standard_subspace_factor_base, KoblitzCurve,
};
use crypto_lib::cryptanalysis::wide_groebner::TwoWordFieldStructure;
use crypto_lib::cryptanalysis::wide_sixsum::{Mono512, RootReduction, System512};
use num_bigint::BigUint;
use serde_json::json;

const PANEL: [usize; 16] = [
    0, 2, 4, 6, 8, 10, 64, 128, 512, 1024, 4096, 8192, 16384, 32768, 49152, 64906,
];

fn x_of(p: &BinaryPoint) -> &F2mElement {
    match p {
        BinaryPoint::Affine { x, .. } => x,
        BinaryPoint::Infinity => panic!("affine required"),
    }
}

fn code_of(x: &F2mElement) -> u64 {
    let words = x.raw_bits();
    assert_eq!(words[0] >> 16, 0);
    assert_eq!(words[1], 0);
    words[0]
}

fn put_element(bits: &mut Mono512, offset: usize, width: usize, x: &F2mElement) {
    let words = x.raw_bits();
    for i in 0..width {
        if words[i / 64] >> (i % 64) & 1 == 1 {
            bits.0[(offset + i) / 64] |= 1u64 << ((offset + i) % 64);
        }
    }
}

fn point_json(p: &BinaryPoint) -> serde_json::Value {
    match p {
        BinaryPoint::Infinity => json!({"infinity":true}),
        BinaryPoint::Affine { x, y } => json!({
            "x":x.to_biguint().to_str_radix(16), "y":y.to_biguint().to_str_radix(16)
        }),
    }
}

fn peak_rss_bytes() -> Option<u64> {
    #[cfg(any(target_os = "macos", target_os = "linux"))]
    {
        let mut usage = std::mem::MaybeUninit::<libc::rusage>::uninit();
        if unsafe { libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) } != 0 {
            return None;
        }
        let usage = unsafe { usage.assume_init() };
        #[cfg(target_os = "macos")]
        return Some(usage.ru_maxrss as u64);
        #[cfg(target_os = "linux")]
        return Some((usage.ru_maxrss as u64) * 1024);
    }
    #[cfg(not(any(target_os = "macos", target_os = "linux")))]
    None
}

fn main() {
    let setup = Instant::now();
    let kc = KoblitzCurve::known_n83_k0().expect("pinned curve");
    assert_eq!(kc.label(), "icv1-f2m83-tm6151469093347-debefd74");
    let source = build_standard_subspace_factor_base(&kc, 16).expect("source base");
    assert_eq!(source.points.len(), 64_907);
    let table = TwoWordFieldStructure::new(83, &kc.curve.irreducible).expect("two-word field");
    let four = BigUint::from(4u32);
    let inv_four = (&kc.subgroup_order * BigUint::from(3u32) + BigUint::from(1u32)) / &four;
    let public_target = BinaryPoint::Affine {
        x: F2mElement::from_hex("355fb5df7a905f16921eb", 83),
        y: F2mElement::from_hex("5900a390f42d290f1bbe", 83),
    };
    assert!(kc.curve.is_on_curve(&public_target));
    assert_eq!(
        kc.mul(&public_target, &kc.subgroup_order),
        BinaryPoint::Infinity
    );
    let source_target = kc.mul(&public_target, &inv_four);
    assert_eq!(kc.mul(&source_target, &four), public_target);
    let ordinary = System512::build(
        &source.subspace_basis,
        x_of(&source_target),
        &kc.curve.b,
        6,
        &table,
    )
    .expect("direct ordinary system");
    assert_eq!(ordinary.n_vars, 428);

    let planted_indices = [0usize, 2, 4, 6, 8, 10];
    let mut sum = BinaryPoint::Infinity;
    let mut assignment = Mono512::default();
    for (i, &index) in planted_indices.iter().enumerate() {
        let point = &source.points[index];
        put_element(&mut assignment, i * 16, 16, x_of(point));
        sum = kc.add(&sum, point);
        if (1..=4).contains(&i) {
            put_element(&mut assignment, 96 + (i - 1) * 83, 83, x_of(&sum));
        }
    }
    let planted = System512::build(&source.subspace_basis, x_of(&sum), &kc.curve.b, 6, &table)
        .expect("direct planted system");
    assert!(planted.all_vanish(&assignment));
    let planted_fixed = planted
        .assign_summand_code(0, 16, code_of(x_of(&source.points[0])))
        .expect("valid first code");
    assert!(planted_fixed.all_vanish(&assignment));
    println!(
        "{}",
        json!({
            "phase":"setup", "curve_id":kc.label(), "source_points":source.points.len(),
            "public_target":point_json(&public_target), "source_preimage":point_json(&source_target),
            "torsion_offset":0, "system_variables":ordinary.n_vars,
            "system_equations":ordinary.equations.len(), "planted_indices":planted_indices,
            "planted_original_vanishes":true, "planted_substituted_vanishes":true,
            "setup_ns":setup.elapsed().as_nanos(), "peak_rss_bytes":peak_rss_bytes(),
            "claim_scope":"prefix_pruning_diagnostic_only"
        })
    );
    std::io::stdout().flush().expect("flush setup receipt");

    let mut seen = HashSet::new();
    for (panel_position, &source_index) in PANEL.iter().enumerate() {
        let x = x_of(&source.points[source_index]);
        let code = code_of(x);
        let duplicate_x = !seen.insert(code);
        let start = Instant::now();
        let fixed = ordinary
            .assign_summand_code(0, 16, code)
            .expect("valid source code");
        let assign_ns = start.elapsed().as_nanos();
        let start = Instant::now();
        let reduction = fixed.root_reduce();
        let reduce_ns = start.elapsed().as_nanos();
        let (status, columns, rank, contradiction, xor_ops, rows) = match reduction {
            RootReduction::ColumnLimit { columns } => {
                ("column_limit", columns, None, None, None, None)
            }
            RootReduction::Reduced {
                columns,
                rank,
                contradiction,
                xor_ops,
                linear_equations,
                ..
            } => {
                let rows: Vec<_> = linear_equations.into_iter().map(|row| {
                    let source_only = row.variables.iter().all(|&v| v < ordinary.summand_bits);
                    json!({"variables":row.variables,"constant":row.constant,"source_only":source_only})
                }).collect();
                (
                    "reduced",
                    columns,
                    Some(rank),
                    Some(contradiction),
                    Some(xor_ops),
                    Some(rows),
                )
            }
        };
        let source_only_rows = rows.as_ref().map(|rs| {
            rs.iter()
                .filter(|r| {
                    r["source_only"] == true
                        && r["variables"].as_array().is_some_and(|v| !v.is_empty())
                })
                .count()
        });
        println!(
            "{}",
            json!({
                "phase":"prefix", "panel_position":panel_position,
                "source_index":source_index, "source_x":x.to_biguint().to_str_radix(16),
                "code":code, "duplicate_x":duplicate_x, "equations":fixed.equations.len(),
                "term_occurrences":fixed.monomial_count(), "status":status,
                "columns":columns,"rank":rank,"contradiction":contradiction,
                "source_only_rows":source_only_rows,"linear_equations":rows,
                "xor_ops":xor_ops,"assign_ns":assign_ns,"reduce_ns":reduce_ns,
                "peak_rss_bytes":peak_rss_bytes(),"claim_scope":"prefix_pruning_diagnostic_only"
            })
        );
        std::io::stdout().flush().expect("flush prefix row");
    }
}
