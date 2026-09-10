//! Export one binary-Koblitz point-decomposition problem from a single
//! source system to WDSat ANF, CryptoMiniSat CNF-XOR, and Magma Boolean F4.
//!
//! The factor base is defined only from field coordinates.  `standard`
//! uses `span(1,z,...,z^(ell-1))`; `ggmp` derives the invariant space from
//! an irreducible factor of `T^n-1`.  Neither path enumerates the target
//! subgroup or attaches discrete-log labels to factor-base points.
//!
//! ```text
//! cargo run --release --example koblitz_pdp_export -- \
//!   31 5 ggmp 20260909 100000 /tmp/pdp-n31-ggmp
//! ```

use crypto_lib::binary_ecc::curve::{point_add, point_neg};
use crypto_lib::binary_ecc::{BinaryCurve, BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::binary_semaev_s4::{weil_descend_s4, S4System};
use crypto_lib::cryptanalysis::koblitz_groebner::{
    build_decomposition_system, DecompositionSystem, FieldStructure,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    factor_x_n_minus_1, find_irreducible_sparse, invariant_subspace_basis, points_with_x,
};
use crypto_lib::cryptanalysis::sat::{SolveResult, SolverStats};
use crypto_lib::cryptanalysis::semaev_sat::{
    encode_boolean_system_with, encode_semaev_s4_with, S4Options, XorEncoding,
};
use num_bigint::BigUint;
use num_traits::{One, Zero};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};
use std::collections::{BTreeMap, HashMap};
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Clone, Debug)]
struct Row {
    /// Zero-based Boolean variable ids in each square-free monomial.
    monomials: Vec<Vec<u32>>,
    constant: bool,
}

fn point_key(point: &BinaryPoint) -> Option<(BigUint, BigUint)> {
    match point {
        BinaryPoint::Infinity => None,
        BinaryPoint::Affine { x, y } => Some((x.to_biguint(), y.to_biguint())),
    }
}

fn rows_from_s4(system: &S4System) -> Vec<Row> {
    let n_x = system.n_x_vars();
    let mut rows = Vec::new();
    for (i, coefficients) in system.correspondence.iter().enumerate() {
        for (degree, sigma) in coefficients.iter().enumerate() {
            let mut monomials: Vec<Vec<u32>> = sigma
                .monomials()
                .filter(|m| !m.is_empty())
                .cloned()
                .collect();
            monomials.push(vec![n_x + system.e_var(i, degree)]);
            rows.push(Row {
                monomials,
                constant: sigma.has_constant(),
            });
        }
    }
    for equation in &system.semaev {
        rows.push(Row {
            monomials: equation
                .monomials()
                .filter(|m| !m.is_empty())
                .map(|m| m.iter().map(|v| n_x + *v).collect())
                .collect(),
            constant: equation.has_constant(),
        });
    }
    rows
}

fn rows_from_chained(system: &DecompositionSystem) -> Vec<Row> {
    system
        .equations
        .iter()
        .map(|equation| Row {
            monomials: equation
                .terms
                .iter()
                .filter(|term| term.mask != 0)
                .map(|term| {
                    (0..system.n_vars)
                        .filter(|v| (term.mask >> v) & 1 == 1)
                        .map(|v| v as u32)
                        .collect()
                })
                .collect(),
            constant: equation.terms.iter().any(|term| term.mask == 0),
        })
        .collect()
}

fn wdsat_anf(n_vars: u32, rows: &[Row]) -> String {
    let active: Vec<&Row> = rows
        .iter()
        .filter(|row| row.constant || !row.monomials.is_empty())
        .collect();
    let mut out = format!("p cnf {n_vars} {}\n", active.len());
    for row in active {
        out.push('x');
        for monomial in &row.monomials {
            if monomial.len() > 1 {
                out.push_str(&format!(" .{}", monomial.len()));
            }
            for variable in monomial {
                out.push_str(&format!(" {}", variable + 1));
            }
        }
        // WDSat's ANF rows have odd parity.  Add T exactly when the
        // polynomial's constant is zero, so the remaining terms equal 0.
        if !row.constant {
            out.push_str(" T");
        }
        out.push_str(" 0\n");
    }
    out
}

fn xor_dimacs(n_vars: u32, rows: &[Row]) -> (String, usize, usize, u32) {
    let mut aux_of: BTreeMap<Vec<u32>, u32> = BTreeMap::new();
    for row in rows {
        for monomial in &row.monomials {
            if monomial.len() >= 2 {
                let mut key = monomial.clone();
                key.sort_unstable();
                key.dedup();
                aux_of.entry(key).or_insert(0);
            }
        }
    }
    let mut next = n_vars + 1;
    for value in aux_of.values_mut() {
        *value = next;
        next += 1;
    }

    let mut clauses: Vec<Vec<i32>> = Vec::new();
    for (variables, auxiliary) in &aux_of {
        let z = *auxiliary as i32;
        let mut forward: Vec<i32> = variables.iter().map(|v| -((*v + 1) as i32)).collect();
        forward.push(z);
        clauses.push(forward);
        for variable in variables {
            clauses.push(vec![(*variable + 1) as i32, -z]);
        }
    }

    let mut xors: Vec<Vec<i32>> = Vec::new();
    for row in rows {
        let mut atoms: Vec<i32> = row
            .monomials
            .iter()
            .map(|monomial| {
                if monomial.len() == 1 {
                    (monomial[0] + 1) as i32
                } else {
                    let mut key = monomial.clone();
                    key.sort_unstable();
                    key.dedup();
                    aux_of[&key] as i32
                }
            })
            .collect();
        if atoms.is_empty() {
            if row.constant {
                clauses.push(Vec::new());
            }
            continue;
        }
        // CryptoMiniSat's extended-DIMACS XOR line has parity 1.
        // Negating one atom flips it to parity 0.
        if !row.constant {
            atoms[0] = -atoms[0];
        }
        xors.push(atoms);
    }

    let constraints = clauses.len() + xors.len();
    let max_var = next - 1;
    let mut out = format!("p cnf {max_var} {constraints}\n");
    for clause in &clauses {
        for literal in clause {
            out.push_str(&format!("{literal} "));
        }
        out.push_str("0\n");
    }
    for xor in &xors {
        out.push('x');
        for literal in xor {
            out.push_str(&format!(" {literal}"));
        }
        out.push_str(" 0\n");
    }
    (out, clauses.len(), xors.len(), max_var)
}

fn magma_script(n_vars: u32, rows: &[Row]) -> String {
    let mut polynomials = Vec::new();
    for row in rows {
        let mut terms: Vec<String> = row
            .monomials
            .iter()
            .map(|monomial| {
                monomial
                    .iter()
                    .map(|v| format!("X[{}]", v + 1))
                    .collect::<Vec<_>>()
                    .join("*")
            })
            .collect();
        if row.constant {
            terms.push("1".to_string());
        }
        if !terms.is_empty() {
            polynomials.push(terms.join(" + "));
        }
    }
    format!(
        "SetNthreads(1);\nR := BooleanPolynomialRing({}, \"grevlex\");\nX := [R.i : i in [1..{}]];\nF := [{}];\nI := ideal<R | F>;\ntime GroebnerBasis(I);\nquit;\n",
        n_vars,
        n_vars,
        polynomials.join(",\n  ")
    )
}

fn point_from_basis(curve: &BinaryCurve, basis: &[F2mElement], rng: &mut StdRng) -> BinaryPoint {
    loop {
        let mut x = F2mElement::zero(curve.m);
        for element in basis {
            if rng.gen::<bool>() {
                x = x.add(element);
            }
        }
        let points = points_with_x(curve, &x);
        if !points.is_empty() {
            return points[rng.gen_range(0..points.len())].clone();
        }
    }
}

fn planted_target(
    curve: &BinaryCurve,
    basis: &[F2mElement],
    seed: u64,
) -> (BinaryPoint, Vec<BinaryPoint>) {
    let mut rng = StdRng::seed_from_u64(seed);
    loop {
        let points: Vec<_> = (0..3)
            .map(|_| point_from_basis(curve, basis, &mut rng))
            .collect();
        let sum = point_add(curve, &point_add(curve, &points[0], &points[1]), &points[2]);
        if sum != BinaryPoint::Infinity {
            return (sum, points);
        }
    }
}

fn direct_mitm(curve: &BinaryCurve, basis: &[F2mElement], target: &BinaryPoint) -> Value {
    if basis.len() > 10 {
        return json!({"status":"not_run","reason":"factor-base materialisation cap","ell_cap":10});
    }
    let start = Instant::now();
    let mut points = Vec::new();
    for mask in 0..(1usize << basis.len()) {
        let mut x = F2mElement::zero(curve.m);
        for (i, element) in basis.iter().enumerate() {
            if (mask >> i) & 1 == 1 {
                x = x.add(element);
            }
        }
        points.extend(points_with_x(curve, &x));
    }
    let setup_ns = start.elapsed().as_nanos();
    let search = Instant::now();
    let mut pairs: HashMap<Option<(BigUint, BigUint)>, (usize, usize)> = HashMap::new();
    let mut additions = 0u64;
    for i in 0..points.len() {
        for j in i..points.len() {
            let sum = point_add(curve, &points[i], &points[j]);
            additions += 1;
            pairs.entry(point_key(&sum)).or_insert((i, j));
        }
    }
    let mut witness = None;
    for (k, point) in points.iter().enumerate() {
        let needed = point_add(curve, target, &point_neg(point));
        additions += 1;
        if let Some(&(i, j)) = pairs.get(&point_key(&needed)) {
            witness = Some((i, j, k));
            break;
        }
    }
    let search_ns = search.elapsed().as_nanos();
    json!({
        "status": if witness.is_some() {"sat"} else {"unsat"},
        "factor_points": points.len(),
        "pair_entries": pairs.len(),
        "group_additions": additions,
        "factor_base_ns": setup_ns,
        "search_ns": search_ns,
        "wall_ns": setup_ns + search_ns,
        "witness_indices": witness.map(|(i,j,k)| vec![i,j,k])
    })
}

fn write(path: &Path, contents: &str) {
    fs::write(path, contents).unwrap_or_else(|error| panic!("write {}: {error}", path.display()));
}

fn stats_json(stats: &SolverStats) -> Value {
    json!({
        "decisions": stats.decisions,
        "conflicts": stats.conflicts,
        "restarts": stats.restarts,
        "propagations": stats.propagations,
        "xor_passes": stats.xor_passes,
        "xor_propagations": stats.xor_propagations,
        "xor_conflicts": stats.xor_conflicts,
        "learnt_clauses": stats.learnt_clauses,
        "max_level": stats.max_level
    })
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert_eq!(
        args.len(),
        7,
        "usage: koblitz_pdp_export <n> <ell> <standard|ggmp> <seed> <conflict-budget> <new-output-dir>"
    );
    let n: u32 = args[1].parse().expect("n");
    let requested_ell: usize = args[2].parse().expect("ell");
    let kind = args[3].as_str();
    let seed: u64 = args[4].parse().expect("seed");
    let conflict_budget: u64 = args[5].parse().expect("conflict budget");
    let output = PathBuf::from(&args[6]);
    assert!(
        !output.exists(),
        "output path must be new: {}",
        output.display()
    );
    fs::create_dir_all(&output).expect("create output directory");

    let total_start = Instant::now();
    let predicate_start = Instant::now();
    let (irreducible, basis, predicate): (_, Vec<F2mElement>, Value) = match kind {
        "standard" => {
            let irr = find_irreducible_sparse(n).expect("sparse irreducible polynomial");
            let basis = (0..requested_ell)
                .map(|i| F2mElement::from_bit_positions(&[i as u32], n))
                .collect();
            (
                irr,
                basis,
                json!({
                    "kind":"polynomial_subspace",
                    "definition":"x = sum(c_i z^i), c_i in F_2, 0 <= i < ell",
                    "enumerates_target_subgroup":false,
                    "uses_discrete_log_labels":false
                }),
            )
        }
        "ggmp" => {
            let factor_index = 0usize;
            let (irr, basis) = invariant_subspace_basis(n, factor_index)
                .expect("GGMP invariant subspace for this degree");
            assert_eq!(
                basis.len(),
                requested_ell,
                "requested ell does not match GGMP factor dimension"
            );
            let factor = factor_x_n_minus_1(n)[factor_index];
            let exponents: Vec<u32> = (0..64).filter(|i| (factor >> i) & 1 == 1).collect();
            (
                irr,
                basis,
                json!({
                    "kind":"ggmp_linearised_kernel",
                    "factor_index":factor_index,
                    "factor_bitmask":factor,
                    "linearised_exponents":exponents,
                    "definition":"kernel of sum X^(2^i) over the set factor exponents",
                    "enumerates_target_subgroup":false,
                    "uses_discrete_log_labels":false
                }),
            )
        }
        _ => panic!("basis must be standard or ggmp"),
    };
    let predicate_ns = predicate_start.elapsed().as_nanos();

    let curve = BinaryCurve {
        m: n,
        irreducible: irreducible.clone(),
        a: F2mElement::one(n),
        b: F2mElement::one(n),
        generator: BinaryPoint::Infinity,
        order: BigUint::zero(),
        cofactor: BigUint::one(),
    };
    let target_start = Instant::now();
    let (target, planted) = planted_target(&curve, &basis, seed);
    let target_ns = target_start.elapsed().as_nanos();
    let target_x = match &target {
        BinaryPoint::Affine { x, .. } => x.clone(),
        BinaryPoint::Infinity => unreachable!(),
    };

    let system_start = Instant::now();
    let (n_vars, rows, representation, chained_system) = if kind == "ggmp" {
        let field = FieldStructure::new(n, &irreducible);
        let system = build_decomposition_system(&basis, &target_x, &F2mElement::one(n), 3, &field)
            .expect("GGMP chained S3 system fits the 64-variable representation");
        (
            system.n_vars as u32,
            rows_from_chained(&system),
            "chained_s3",
            Some(system),
        )
    } else {
        let system = weil_descend_s4(
            n,
            requested_ell as u32,
            &irreducible,
            &F2mElement::one(n),
            &target_x,
        );
        (
            system.n_x_vars() + system.n_e_vars(),
            rows_from_s4(&system),
            "symmetrised_s4",
            None,
        )
    };
    let system_ns = system_start.elapsed().as_nanos();

    let anf = wdsat_anf(n_vars, &rows);
    let (cms, cnf_clauses, xor_rows, cms_vars) = xor_dimacs(n_vars, &rows);
    let magma = magma_script(n_vars, &rows);
    write(&output.join("instance.anf"), &anf);
    write(&output.join("instance.xor.cnf"), &cms);
    write(&output.join("instance.magma"), &magma);

    let native_start = Instant::now();
    let (native_result, native_stats, native_vars, native_clauses, model_valid) =
        if let Some(system) = chained_system {
            let mut encoding = encode_boolean_system_with(
                system.n_vars,
                &system.equations,
                &[],
                XorEncoding::Native,
            );
            encoding.solver.conflict_budget = conflict_budget;
            let initial_clauses = encoding.solver.n_clauses();
            let result = encoding.solver.solve();
            let valid = (result == SolveResult::Sat).then(|| {
                let assignment = encoding.solver.model();
                system.equations.iter().all(|equation| {
                    equation.terms.iter().fold(false, |parity, term| {
                        let monomial = if term.mask == 0 {
                            true
                        } else {
                            (0..system.n_vars)
                                .filter(|i| (term.mask >> i) & 1 == 1)
                                .all(|i| assignment[i])
                        };
                        parity ^ monomial
                    }) == false
                })
            });
            (
                format!("{result:?}").to_lowercase(),
                stats_json(&encoding.solver.stats),
                encoding.solver.n_vars(),
                initial_clauses,
                valid,
            )
        } else {
            let mut encoding = encode_semaev_s4_with(
                n,
                requested_ell as u32,
                &irreducible,
                &F2mElement::one(n),
                &target_x,
                S4Options {
                    encoding: XorEncoding::Native,
                    // The external ANF/CNF-XOR/Magma files below carry
                    // exactly the source equations and no ordering clauses.
                    // Keep this arm byte-semantically matched to them.
                    break_symmetry: false,
                },
            );
            encoding.solver.conflict_budget = conflict_budget;
            let initial_clauses = encoding.solver.n_clauses();
            let result = encoding.solver.solve();
            let valid = (result == SolveResult::Sat).then(|| {
                let decoded = encoding.decode();
                let check = crypto_lib::cryptanalysis::binary_semaev_s4::elementary_symmetric_3(
                    &decoded[0],
                    &decoded[1],
                    &decoded[2],
                    &irreducible,
                );
                crypto_lib::cryptanalysis::binary_semaev_s4::symmetrised_s4_eval(
                    &check.0,
                    &check.1,
                    &check.2,
                    &target_x,
                    &irreducible,
                )
                .is_zero()
            });
            (
                format!("{result:?}").to_lowercase(),
                stats_json(&encoding.solver.stats),
                encoding.solver.n_vars(),
                initial_clauses,
                valid,
            )
        };
    let native_ns = native_start.elapsed().as_nanos();
    assert_ne!(
        model_valid,
        Some(false),
        "native SAT model failed the source ANF system"
    );

    let mitm = direct_mitm(&curve, &basis, &target);
    let basis_bits: Vec<String> = basis.iter().map(|x| x.to_biguint().to_string()).collect();
    let planted_points: Vec<Value> = planted
        .iter()
        .map(|point| match point {
            BinaryPoint::Affine { x, y } => json!({
                "x":x.to_biguint().to_string(),
                "y":y.to_biguint().to_string()
            }),
            BinaryPoint::Infinity => Value::Null,
        })
        .collect();
    let target_json = match &target {
        BinaryPoint::Affine { x, y } => json!({
            "x":x.to_biguint().to_string(),
            "y":y.to_biguint().to_string()
        }),
        BinaryPoint::Infinity => Value::Null,
    };
    let manifest = json!({
        "kind":"binary_koblitz_pdp_cross_solver_instance",
        "scope":"standalone planted point-decomposition problem; not a completed index-calculus attack",
        "n":n,
        "ell":basis.len(),
        "m":3,
        "seed":seed,
        "curve":"y^2 + xy = x^3 + x^2 + 1",
        "irreducible_low_terms":irreducible.low_terms,
        "factor_base_predicate":predicate,
        "factor_base_basis_bitmasks":basis_bits,
        "target":target_json,
        "planted_points":planted_points,
        "representation":representation,
        "source_variables":n_vars,
        "source_equations":rows.len(),
        "source_max_degree":rows.iter().flat_map(|row| &row.monomials).map(Vec::len).max().unwrap_or(0),
        "source_max_monomials_per_equation":rows.iter().map(|row| row.monomials.len()).max().unwrap_or(0),
        "exports":{
            "wdsat_anf":{"path":"instance.anf","bytes":anf.len(),"blake3":blake3::hash(anf.as_bytes()).to_hex().to_string()},
            "cryptominisat_xor_dimacs":{"path":"instance.xor.cnf","variables":cms_vars,"cnf_clauses":cnf_clauses,"xor_rows":xor_rows,"bytes":cms.len(),"blake3":blake3::hash(cms.as_bytes()).to_hex().to_string()},
            "magma_boolean_f4":{"path":"instance.magma","bytes":magma.len(),"blake3":blake3::hash(magma.as_bytes()).to_hex().to_string()}
        },
        "timing_ns":{
            "factor_base_predicate_construction":predicate_ns,
            "planted_target_construction":target_ns,
            "source_system_construction":system_ns,
            "native_encoding_and_solve":native_ns,
            "whole_process_internal":total_start.elapsed().as_nanos()
        },
        "native_sat":{
            "result":native_result,
            "conflict_budget":conflict_budget,
            "variables":native_vars,
            "initial_cnf_clauses":native_clauses,
            "source_model_valid":model_valid,
            "stats":native_stats
        },
        "direct_meet_in_the_middle":mitm,
        "interpretation":"SAT is a witness only after source-system validation; Unknown is inconclusive and is never counted as UNSAT"
    });
    let manifest_text = serde_json::to_string_pretty(&manifest).expect("manifest json") + "\n";
    write(&output.join("manifest.json"), &manifest_text);
    println!("{manifest_text}");
}
