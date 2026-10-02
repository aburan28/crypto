//! Exact search over nonlinear label gauges of a three-bit factor-base cube.

use crypto_lib::binary_ecc::F2mElement;
use crypto_lib::cryptanalysis::koblitz_groebner::{
    build_decomposition_system, solving_profile_sparse, system_degree, FieldStructure,
    SolvingProfile,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::find_irreducible_sparse;
use crypto_lib::cryptanalysis::pq_groebner_f2::{F2BoolMono, F2BoolPoly};
use num_bigint::BigUint;
use serde::Serialize;
use serde_json::json;
use std::cmp::Ordering;
use std::collections::{BTreeSet, HashSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

const N: u32 = 7;
const ELL: usize = 3;
const M: usize = 3;
const LEAF_VARS: usize = ELL * M;

#[derive(Clone, Copy)]
struct FrozenCase {
    split: &'static str,
    draw: u32,
    basis: [u64; ELL],
    target: u64,
}

const CASES: [FrozenCase; 4] = [
    FrozenCase {
        split: "discovery",
        draw: 3,
        basis: [22, 30, 88],
        target: 9,
    },
    FrozenCase {
        split: "discovery",
        draw: 4,
        basis: [35, 87, 26],
        target: 24,
    },
    FrozenCase {
        split: "holdout",
        draw: 9,
        basis: [120, 70, 119],
        target: 34,
    },
    FrozenCase {
        split: "holdout",
        draw: 11,
        basis: [12, 69, 104],
        target: 7,
    },
];

#[derive(Clone, Serialize)]
struct Profile {
    degree: u32,
    rows: usize,
    cols: usize,
    rank: usize,
    refuted: bool,
    vars_determined: usize,
    vars_occurring: usize,
    resolves: bool,
}

impl From<SolvingProfile> for Profile {
    fn from(p: SolvingProfile) -> Self {
        Self {
            degree: p.degree,
            rows: p.rows,
            cols: p.cols,
            rank: p.rank,
            refuted: p.refuted,
            vars_determined: p.vars_determined,
            vars_occurring: p.vars_occurring,
            resolves: p.resolves(),
        }
    }
}

#[derive(Clone, Serialize)]
struct CaseResult {
    split: String,
    draw: u32,
    input_degree: u32,
    input_terms: usize,
    resolution_degree_through_6: Option<u32>,
    degree_5: Option<Profile>,
    degree_6: Option<Profile>,
    exhaustive_equivalence: bool,
    exhaustive_solution_count: usize,
    elapsed_ms: u128,
}

#[derive(Clone, Serialize)]
struct CandidateResult {
    class_index: usize,
    permutation: [u8; 8],
    permutation_code: u32,
    anf_degree: u32,
    affine: bool,
    cases: Vec<CaseResult>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct SelectionKey {
    unresolved_through_5: usize,
    max_resolution_through_6: usize,
    degree_5_cols: usize,
    degree_5_rows: usize,
    input_terms: usize,
    permutation_code: u32,
}

impl Ord for SelectionKey {
    fn cmp(&self, other: &Self) -> Ordering {
        (
            self.unresolved_through_5,
            self.max_resolution_through_6,
            self.degree_5_cols,
            self.degree_5_rows,
            self.input_terms,
            self.permutation_code,
        )
            .cmp(&(
                other.unresolved_through_5,
                other.max_resolution_through_6,
                other.degree_5_cols,
                other.degree_5_rows,
                other.input_terms,
                other.permutation_code,
            ))
    }
}

impl PartialOrd for SelectionKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

fn encode_perm(p: &[u8; 8]) -> u32 {
    p.iter()
        .enumerate()
        .fold(0u32, |word, (i, &x)| word | (u32::from(x) << (3 * i)))
}

fn next_permutation(p: &mut [u8; 8]) -> bool {
    let Some(i) = (0..p.len() - 1).rfind(|&i| p[i] < p[i + 1]) else {
        return false;
    };
    let j = (i + 1..p.len()).rfind(|&j| p[i] < p[j]).unwrap();
    p.swap(i, j);
    p[i + 1..].reverse();
    true
}

fn linear_image(rows: [u8; 3], x: u8) -> u8 {
    rows.iter().enumerate().fold(0, |y, (i, row)| {
        y | (((row & x).count_ones() as u8 & 1) << i)
    })
}

fn affine_group() -> Vec<[u8; 8]> {
    let mut group = BTreeSet::new();
    for packed in 0u16..512 {
        let rows = [
            (packed & 7) as u8,
            ((packed >> 3) & 7) as u8,
            ((packed >> 6) & 7) as u8,
        ];
        let linear: Vec<_> = (0..8).map(|x| linear_image(rows, x)).collect();
        if linear.iter().copied().collect::<BTreeSet<_>>().len() != 8 {
            continue;
        }
        for offset in 0..8 {
            let mut p = [0u8; 8];
            for x in 0..8 {
                p[x] = linear[x] ^ offset;
            }
            group.insert(p);
        }
    }
    group.into_iter().collect()
}

fn right_coset_representatives(affine: &[[u8; 8]]) -> Vec<[u8; 8]> {
    let mut covered = HashSet::new();
    let mut reps = Vec::new();
    let mut p = [0, 1, 2, 3, 4, 5, 6, 7];
    loop {
        if !covered.contains(&encode_perm(&p)) {
            reps.push(p);
            for a in affine {
                let mut q = [0u8; 8];
                for x in 0..8 {
                    q[x] = p[a[x] as usize];
                }
                covered.insert(encode_perm(&q));
            }
        }
        if !next_permutation(&mut p) {
            break;
        }
    }
    assert_eq!(covered.len(), 40_320, "affine right cosets must cover S8");
    reps
}

/// ANF monomial masks for the three output bits of `p(new_label)`.
fn anf_forms(p: &[u8; 8]) -> [Vec<u8>; 3] {
    std::array::from_fn(|bit| {
        let mut truth = [0u8; 8];
        for x in 0..8 {
            truth[x] = (p[x] >> bit) & 1;
        }
        for variable in 0..3 {
            for mask in 0..8 {
                if mask & (1 << variable) != 0 {
                    truth[mask] ^= truth[mask ^ (1 << variable)];
                }
            }
        }
        (0..8).filter(|&mask| truth[mask] == 1).map(|x| x as u8).collect()
    })
}

fn anf_degree(p: &[u8; 8]) -> u32 {
    anf_forms(p)
        .iter()
        .flatten()
        .map(|m| m.count_ones())
        .max()
        .unwrap_or(0)
}

fn toggle(set: &mut BTreeSet<u64>, monomial: u64) {
    if !set.insert(monomial) {
        set.remove(&monomial);
    }
}

fn multiply_masks(a: &BTreeSet<u64>, b: &[u64]) -> BTreeSet<u64> {
    let mut out = BTreeSet::new();
    for &x in a {
        for &y in b {
            toggle(&mut out, x | y);
        }
    }
    out
}

/// Simultaneously replace each old leaf bit by the ANF of its new label.
fn transform_poly(poly: &F2BoolPoly, p: &[u8; 8]) -> F2BoolPoly {
    let local = anf_forms(p);
    let forms: Vec<Vec<u64>> = (0..LEAF_VARS)
        .map(|v| {
            let block = v / ELL;
            let coordinate = v % ELL;
            local[coordinate]
                .iter()
                .map(|&mask| u64::from(mask) << (block * ELL))
                .collect()
        })
        .collect();
    let leaf_mask = (1u64 << LEAF_VARS) - 1;
    let mut out = BTreeSet::new();
    for term in &poly.terms {
        let mut expansion = BTreeSet::from([term.mask & !leaf_mask]);
        let mut leaves = term.mask & leaf_mask;
        while leaves != 0 {
            let v = leaves.trailing_zeros() as usize;
            leaves &= leaves - 1;
            expansion = multiply_masks(&expansion, &forms[v]);
        }
        for monomial in expansion {
            toggle(&mut out, monomial);
        }
    }
    F2BoolPoly::from_monos(
        out.into_iter().map(F2BoolMono::from_mask).collect(),
        poly.n_vars,
    )
}

fn transform_system(polys: &[F2BoolPoly], p: &[u8; 8]) -> Vec<F2BoolPoly> {
    polys.iter().map(|poly| transform_poly(poly, p)).collect()
}

fn eval(poly: &F2BoolPoly, assignment: u64) -> bool {
    poly.terms
        .iter()
        .filter(|term| term.mask & !assignment == 0)
        .count()
        % 2
        == 1
}

fn mapped_assignment(new_assignment: u64, p: &[u8; 8]) -> u64 {
    let mut old = new_assignment & !((1u64 << LEAF_VARS) - 1);
    for block in 0..M {
        let code = ((new_assignment >> (block * ELL)) & 7) as usize;
        old |= u64::from(p[code]) << (block * ELL);
    }
    old
}

fn exhaustive_equivalence(
    original: &[F2BoolPoly],
    transformed: &[F2BoolPoly],
    n_vars: usize,
    p: &[u8; 8],
) -> (bool, usize) {
    let mut solutions = 0;
    for assignment in 0..(1u64 << n_vars) {
        let mapped = mapped_assignment(assignment, p);
        let mut all_zero = true;
        for (a, b) in original.iter().zip(transformed) {
            let expected = eval(a, mapped);
            let got = eval(b, assignment);
            if expected != got {
                return (false, solutions);
            }
            all_zero &= !got;
        }
        solutions += usize::from(all_zero);
    }
    (true, solutions)
}

fn inverse_permutation(p: &[u8; 8]) -> [u8; 8] {
    let mut inverse = [0u8; 8];
    for (x, &y) in p.iter().enumerate() {
        inverse[y as usize] = x as u8;
    }
    inverse
}

fn verify_inverse(p: &[u8; 8]) {
    let inverse = inverse_permutation(p);
    for x in 0..8 {
        assert_eq!(inverse[p[x] as usize], x as u8);
        assert_eq!(p[inverse[x] as usize], x as u8);
    }
}

fn element(word: u64) -> F2mElement {
    F2mElement::from_biguint(&BigUint::from(word), N)
}

fn build_case(case: FrozenCase, field: &FieldStructure) -> (Vec<F2BoolPoly>, usize) {
    let basis: Vec<_> = case.basis.into_iter().map(element).collect();
    let system = build_decomposition_system(
        &basis,
        &element(case.target),
        &F2mElement::one(N),
        M,
        field,
    )
    .expect("frozen system fits the 64-variable engine");
    assert_eq!(system.n_vars, 16);
    assert_eq!(system.equations.len(), 14);
    (system.equations, system.n_vars)
}

fn measure_case(
    case: FrozenCase,
    original: &[F2BoolPoly],
    n_vars: usize,
    p: &[u8; 8],
) -> CaseResult {
    let started = Instant::now();
    let transformed = transform_system(original, p);
    let input_degree = system_degree(&transformed);
    let input_terms = transformed.iter().map(|poly| poly.terms.len()).sum();
    let mut degree_5 = None;
    let mut degree_6 = None;
    let mut resolution = None;
    for degree in input_degree.max(1)..=6 {
        let profile = solving_profile_sparse(&transformed, n_vars, degree).map(Profile::from);
        if let Some(profile) = &profile {
            if resolution.is_none() && profile.resolves {
                resolution = Some(degree);
            }
        }
        if degree == 5 {
            degree_5 = profile.clone();
        }
        if degree == 6 {
            degree_6 = profile;
        }
        if resolution.is_some() && degree < 5 {
            // Still build degree 5 because it is a frozen selection metric.
            continue;
        }
    }
    let (equivalent, solutions) = exhaustive_equivalence(original, &transformed, n_vars, p);
    CaseResult {
        split: case.split.to_owned(),
        draw: case.draw,
        input_degree,
        input_terms,
        resolution_degree_through_6: resolution,
        degree_5,
        degree_6,
        exhaustive_equivalence: equivalent,
        exhaustive_solution_count: solutions,
        elapsed_ms: started.elapsed().as_millis(),
    }
}

fn selection_key(result: &CandidateResult) -> SelectionKey {
    let unresolved_through_5 = result
        .cases
        .iter()
        .filter(|case| !matches!(case.resolution_degree_through_6, Some(d) if d <= 5))
        .count();
    let max_resolution_through_6 = result
        .cases
        .iter()
        .map(|case| case.resolution_degree_through_6.unwrap_or(99) as usize)
        .max()
        .unwrap_or(99);
    let degree_5_cols = result
        .cases
        .iter()
        .map(|case| case.degree_5.as_ref().map(|p| p.cols).unwrap_or(1usize << 50))
        .sum();
    let degree_5_rows = result
        .cases
        .iter()
        .map(|case| case.degree_5.as_ref().map(|p| p.rows).unwrap_or(1usize << 50))
        .sum();
    SelectionKey {
        unresolved_through_5,
        max_resolution_through_6,
        degree_5_cols,
        degree_5_rows,
        input_terms: result.cases.iter().map(|case| case.input_terms).sum(),
        permutation_code: result.permutation_code,
    }
}

fn measure_candidate(
    class_index: usize,
    p: [u8; 8],
    cases: &[(FrozenCase, Vec<F2BoolPoly>, usize)],
) -> CandidateResult {
    verify_inverse(&p);
    let anf_degree = anf_degree(&p);
    let measured = cases
        .iter()
        .map(|(case, system, n_vars)| measure_case(*case, system, *n_vars, &p))
        .collect();
    CandidateResult {
        class_index,
        permutation: p,
        permutation_code: encode_perm(&p),
        anf_degree,
        affine: anf_degree <= 1,
        cases: measured,
    }
}

fn write_json(path: &Path, value: &impl Serialize) {
    let bytes = serde_json::to_vec_pretty(value).expect("serialize result");
    fs::write(path, bytes).expect("write result");
}

fn output_directory() -> PathBuf {
    let args: Vec<_> = std::env::args().collect();
    assert_eq!(args.len(), 3, "usage: nonlinear-factor-labels --output DIR");
    assert_eq!(args[1], "--output");
    PathBuf::from(&args[2])
}

fn main() {
    let output = output_directory();
    assert!(!output.exists(), "refusing to overwrite {}", output.display());
    fs::create_dir_all(&output).expect("create output directory");
    let started = Instant::now();

    let affine = affine_group();
    assert_eq!(affine.len(), 1_344, "|AGL(3,2)|");
    let representatives = right_coset_representatives(&affine);
    assert_eq!(representatives.len(), 30, "S8 / AGL(3,2) right cosets");
    assert_eq!(representatives[0], [0, 1, 2, 3, 4, 5, 6, 7]);

    let modulus = find_irreducible_sparse(N).expect("degree-7 irreducible modulus");
    let field = FieldStructure::new(N, &modulus);
    let frozen: Vec<_> = CASES
        .iter()
        .map(|case| {
            let (system, n_vars) = build_case(*case, &field);
            (*case, system, n_vars)
        })
        .collect();
    let discovery = &frozen[..2];
    let holdout = &frozen[2..];

    let mut discovery_results = Vec::new();
    for (class_index, &p) in representatives.iter().enumerate() {
        let result = measure_candidate(class_index, p, discovery);
        assert!(result.cases.iter().all(|case| {
            case.exhaustive_equivalence && case.exhaustive_solution_count == 0
        }));
        eprintln!(
            "class {class_index:02}/29 anf={} key={:?} elapsed={}ms",
            result.anf_degree,
            selection_key(&result),
            result.cases.iter().map(|case| case.elapsed_ms).sum::<u128>()
        );
        discovery_results.push(result);
    }
    write_json(&output.join("discovery.json"), &discovery_results);

    let baseline = discovery_results
        .iter()
        .find(|result| result.affine)
        .expect("affine identity class");
    assert_eq!(baseline.class_index, 0);
    assert!(baseline
        .cases
        .iter()
        .all(|case| case.resolution_degree_through_6 == Some(6)));
    let winner = discovery_results
        .iter()
        .filter(|result| !result.affine)
        .min_by_key(|result| selection_key(result))
        .expect("nonlinear class winner")
        .clone();
    let winner_record = json!({
        "frozen_before_holdout": true,
        "class_index": winner.class_index,
        "permutation": winner.permutation,
        "permutation_code": winner.permutation_code,
        "anf_degree": winner.anf_degree,
        "selection_key": format!("{:?}", selection_key(&winner)),
        "discovery_cases": winner.cases,
    });
    write_json(&output.join("winner.json"), &winner_record);

    let baseline_holdout = measure_candidate(0, representatives[0], holdout);
    let winner_holdout = measure_candidate(winner.class_index, winner.permutation, holdout);
    assert!(baseline_holdout.cases.iter().all(|case| {
        case.exhaustive_equivalence
            && case.exhaustive_solution_count == 0
            && case.resolution_degree_through_6 == Some(6)
    }));
    assert!(winner_holdout
        .cases
        .iter()
        .all(|case| case.exhaustive_equivalence && case.exhaustive_solution_count == 0));

    let primary_success = winner_holdout
        .cases
        .iter()
        .all(|case| matches!(case.resolution_degree_through_6, Some(d) if d <= 5));
    let secondary_success = !primary_success
        && winner_holdout.cases.iter().zip(&baseline_holdout.cases).all(|(winner, baseline)| {
            winner.resolution_degree_through_6 == baseline.resolution_degree_through_6
                && winner
                    .degree_6
                    .as_ref()
                    .zip(baseline.degree_6.as_ref())
                    .is_some_and(|(w, b)| w.cols < b.cols)
        });
    let result = json!({
        "schema_version": 1,
        "scope": "GF(2^7) bounded algebra-stage diagnostic; not ECC2K-130 evidence",
        "search_space": {
            "bijections": 40320,
            "affine_group": affine.len(),
            "right_cosets": representatives.len()
        },
        "field_modulus": {
            "degree": modulus.degree,
            "low_terms": modulus.low_terms,
        },
        "baseline_discovery": baseline,
        "winner_discovery": winner,
        "baseline_holdout": baseline_holdout,
        "winner_holdout": winner_holdout,
        "primary_success": primary_success,
        "secondary_success": secondary_success,
        "elapsed_ms": started.elapsed().as_millis(),
        "full_ic_cost": null,
        "rho_ratio": null,
        "m83_gate": null,
        "gf2_131_result": null
    });
    write_json(&output.join("results.json"), &result);
    println!("{}", serde_json::to_string_pretty(&result).unwrap());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn affine_group_and_cosets_have_registered_sizes() {
        let affine = affine_group();
        assert_eq!(affine.len(), 1_344);
        let reps = right_coset_representatives(&affine);
        assert_eq!(reps.len(), 30);
        assert!(reps[0].iter().copied().eq(0..8));
        assert_eq!(reps.iter().filter(|p| anf_degree(p) <= 1).count(), 1);
    }

    #[test]
    fn anf_reconstructs_every_permutation_representative() {
        let reps = right_coset_representatives(&affine_group());
        for p in reps {
            let forms = anf_forms(&p);
            for x in 0u8..8 {
                let mut y = 0u8;
                for (bit, form) in forms.iter().enumerate() {
                    let value = form
                        .iter()
                        .filter(|&&monomial| monomial & !x == 0)
                        .count()
                        % 2;
                    y |= (value as u8) << bit;
                }
                assert_eq!(y, p[x as usize]);
            }
        }
    }

    #[test]
    fn identity_transform_is_exact() {
        let modulus = find_irreducible_sparse(N).unwrap();
        let field = FieldStructure::new(N, &modulus);
        let (system, _) = build_case(CASES[0], &field);
        assert_eq!(transform_system(&system, &[0, 1, 2, 3, 4, 5, 6, 7]), system);
    }

    #[test]
    fn nonlinear_transform_commutes_with_all_assignments() {
        let modulus = find_irreducible_sparse(N).unwrap();
        let field = FieldStructure::new(N, &modulus);
        let (system, n_vars) = build_case(CASES[0], &field);
        let reps = right_coset_representatives(&affine_group());
        let p = reps[1];
        let transformed = transform_system(&system, &p);
        assert_eq!(exhaustive_equivalence(&system, &transformed, n_vars, &p), (true, 0));
    }
}
