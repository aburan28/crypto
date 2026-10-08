//! Five-summand dimension-18 gate: research/f6_n83_fivesum_d18_20261005.
use std::io::Write;
use std::time::Instant;

use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_standard_subspace_factor_base, cofactor_project_factor_base, KoblitzCurve,
};
use crypto_lib::cryptanalysis::wide_groebner::TwoWordFieldStructure;
use crypto_lib::cryptanalysis::wide_sixsum::{Mono512, RootReduction, System512};
use num_bigint::BigUint;
use serde_json::json;

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

fn x_of(point: &BinaryPoint) -> &F2mElement {
    match point {
        BinaryPoint::Affine { x, .. } => x,
        BinaryPoint::Infinity => panic!("affine point required"),
    }
}

fn set_element(bits: &mut Mono512, offset: usize, width: usize, element: &F2mElement) {
    let words = element.raw_bits();
    for i in 0..width {
        if words[i / 64] >> (i % 64) & 1 == 1 {
            bits.0[(offset + i) / 64] |= 1u64 << ((offset + i) % 64);
        }
    }
}

fn point_json(point: &BinaryPoint) -> serde_json::Value {
    match point {
        BinaryPoint::Infinity => json!({"infinity":true}),
        BinaryPoint::Affine { x, y } => json!({
            "x":x.to_biguint().to_str_radix(16), "y":y.to_biguint().to_str_radix(16)
        }),
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert!(
        args.len() == 2 || args.len() == 3,
        "usage: probe planted | ordinary OFFSET(0..3)"
    );
    let start = Instant::now();
    let kc = KoblitzCurve::known_n83_k0().expect("pinned n83 K0");
    assert_eq!(kc.label(), "icv1-f2m83-tm6151469093347-debefd74");
    let source = build_standard_subspace_factor_base(&kc, 18).expect("standard dimension-18 base");
    let projected = cofactor_project_factor_base(&kc, &source).expect("usable subgroup base");
    let usable_points = projected.points.len();
    let signed_columns = projected.signed_orbits.len();
    drop(projected);
    let table = TwoWordFieldStructure::new(83, &kc.curve.irreducible).expect("field table");
    let four = BigUint::from(4u32);
    let inv_four = (&kc.subgroup_order * BigUint::from(3u32) + BigUint::from(1u32)) / &four;
    let one = F2mElement::one(83);
    let zero = F2mElement::zero(83);
    let torsion = [
        BinaryPoint::Infinity,
        BinaryPoint::Affine {
            x: zero.clone(),
            y: one.clone(),
        },
        BinaryPoint::Affine {
            x: one.clone(),
            y: zero,
        },
        BinaryPoint::Affine {
            x: one.clone(),
            y: one,
        },
    ];
    for point in &torsion {
        assert!(kc.curve.is_on_curve(point));
        assert_eq!(kc.mul(point, &four), BinaryPoint::Infinity);
    }
    let public_target = BinaryPoint::Affine {
        x: F2mElement::from_hex("355fb5df7a905f16921eb", 83),
        y: F2mElement::from_hex("5900a390f42d290f1bbe", 83),
    };
    assert!(kc.curve.is_on_curve(&public_target));
    assert_eq!(
        kc.mul(&public_target, &kc.subgroup_order),
        BinaryPoint::Infinity
    );
    let setup_ns = start.elapsed().as_nanos();

    let (target, planted, offset) = match args[1].as_str() {
        "planted" => {
            let chosen = [0usize, 2, 4, 6, 8];
            let sum = chosen.iter().fold(BinaryPoint::Infinity, |acc, &i| {
                kc.add(&acc, &source.points[i])
            });
            assert_ne!(sum, BinaryPoint::Infinity);
            (sum, Some(chosen), None)
        }
        "ordinary" => {
            let i: usize = args
                .get(2)
                .expect("offset required")
                .parse()
                .expect("integer offset");
            assert!(i < 4);
            let subgroup_preimage = kc.mul(&public_target, &inv_four);
            assert_eq!(kc.mul(&subgroup_preimage, &four), public_target);
            (kc.add(&subgroup_preimage, &torsion[i]), None, Some(i))
        }
        _ => panic!("unknown mode"),
    };

    let build = Instant::now();
    let system = System512::build(
        &source.subspace_basis,
        x_of(&target),
        &kc.curve.b,
        5,
        &table,
    )
    .expect("339-variable system admitted");
    let build_ns = build.elapsed().as_nanos();
    assert_eq!(system.n_vars, 339);
    let mut planted_control = None;
    if let Some(chosen) = planted {
        let mut assignment = Mono512::default();
        let mut prefix = BinaryPoint::Infinity;
        for (i, &index) in chosen.iter().enumerate() {
            let point = &source.points[index];
            set_element(&mut assignment, i * 18, 18, x_of(point));
            prefix = kc.add(&prefix, point);
            if (1..=3).contains(&i) {
                set_element(&mut assignment, 90 + (i - 1) * 83, 83, x_of(&prefix));
            }
        }
        assert_eq!(prefix, target);
        assert!(system.all_vanish(&assignment));
        let projected = kc.mul(&target, &four);
        let subgroup_preimage = kc.mul(&projected, &inv_four);
        assert_eq!(kc.mul(&subgroup_preimage, &four), projected);
        let matching_offsets: Vec<_> = torsion
            .iter()
            .enumerate()
            .filter(|(_, t)| kc.add(&subgroup_preimage, t) == target)
            .map(|(i, _)| i)
            .collect();
        assert_eq!(matching_offsets.len(), 1);
        planted_control = Some(json!({
            "indices":chosen, "source_target":point_json(&target),
            "projected_target":point_json(&projected), "all_equations_zero":true,
            "group_sum_replayed":true, "torsion_offset":matching_offsets[0]
        }));
    }
    println!(
        "{}",
        json!({
            "phase":"constructed", "mode":args[1], "offset":offset,
            "curve_id":kc.label(), "source_points":source.points.len(),
            "usable_points":usable_points, "signed_columns":signed_columns,
            "target":point_json(&target), "n_vars":system.n_vars,
            "equations":system.equations.len(), "term_occurrences":system.monomial_count(),
            "max_degree":system.max_degree(), "setup_ns":setup_ns, "build_ns":build_ns,
            "peak_rss_bytes":peak_rss_bytes(), "planted_control":planted_control,
            "claim_scope":"algebraic_feasibility_only"
        })
    );
    std::io::stdout()
        .flush()
        .expect("flush construction receipt");
    if offset.is_some() {
        let reduce = Instant::now();
        let result = system.root_reduce();
        let (status, columns, rank, linear, linear_equations, contradiction, xor_ops) = match result
        {
            RootReduction::ColumnLimit { columns } => {
                ("column_limit", columns, None, None, None, None, None)
            }
            RootReduction::Reduced {
                columns,
                rank,
                linear,
                linear_equations,
                contradiction,
                xor_ops,
            } => (
                "reduced",
                columns,
                Some(rank),
                Some(linear),
                Some(
                    linear_equations
                        .into_iter()
                        .map(|equation| {
                            let summand_variables = equation
                                .variables
                                .iter()
                                .filter(|&&v| v < system.summand_bits)
                                .count();
                            json!({
                                "variables":equation.variables,
                                "constant":equation.constant,
                                "summand_variables":summand_variables,
                                "source_only":summand_variables == equation.variables.len()
                            })
                        })
                        .collect::<Vec<_>>(),
                ),
                Some(contradiction),
                Some(xor_ops),
            ),
        };
        println!(
            "{}",
            json!({
                "phase":"root_reduction", "mode":args[1], "offset":offset,
                "status":status, "columns":columns, "rank":rank,
                "linear":linear, "linear_equations":linear_equations,
                "contradiction":contradiction, "xor_ops":xor_ops,
                "reduce_ns":reduce.elapsed().as_nanos(), "peak_rss_bytes":peak_rss_bytes(),
                "claim_scope":"algebraic_feasibility_only"
            })
        );
    }
}
