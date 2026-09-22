//! Exact reference checks and opt-in, paired kernel measurements.
use super::*;
use serde_json::json;
use std::{hint::black_box, time::Instant};

// Frozen verbatim from 466e8b651002ad11e5de996619f9aaaedb5ac551,
// except for the name. Keep the old indexing and full-width XOR loop.
fn reference(matrix: &mut [Vec<u64>], n_cols: usize, word_ops: &mut u64) -> usize {
    let words = n_cols.div_ceil(64);
    let mut pivot_row = 0usize;
    for c in 0..n_cols {
        let (w, bit) = (c / 64, 1u64 << (c % 64));
        let piv = (pivot_row..matrix.len()).find(|&r| matrix[r][w] & bit != 0);
        let piv = match piv {
            Some(p) => p,
            None => continue,
        };
        matrix.swap(pivot_row, piv);
        for r in 0..matrix.len() {
            if r != pivot_row && matrix[r][w] & bit != 0 {
                for k in 0..words {
                    matrix[r][k] ^= matrix[pivot_row][k];
                }
                *word_ops += words as u64;
            }
        }
        pivot_row += 1;
        if pivot_row == matrix.len() {
            break;
        }
    }
    pivot_row
}

fn next(state: &mut u64) -> u64 {
    // SplitMix64 output mixing avoids the rank-64 cap of a linear xorshift
    // stream when successive output words are arranged as matrix rows.
    *state = state.wrapping_add(0x9e3779b97f4a7c15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
    z ^ (z >> 31)
}

fn synthetic(rows: usize, cols: usize, kind: &str, mut seed: u64) -> Vec<Vec<u64>> {
    let mut out = vec![vec![0; cols.div_ceil(64)]; rows];
    for row in &mut out {
        for word in row.iter_mut() {
            *word = next(&mut seed);
            if kind == "sparse" {
                *word &= next(&mut seed) & next(&mut seed) & next(&mut seed);
            }
        }
        if cols % 64 != 0 {
            *row.last_mut().unwrap() &= (1u64 << (cols % 64)) - 1;
        }
    }
    if kind == "deficient" {
        for r in rows / 2..rows {
            out[r] = out[r - rows / 2].clone();
        }
    }
    out
}

fn check(matrix: &[Vec<u64>], cols: usize) {
    let mut a = matrix.to_vec();
    let mut b = matrix.to_vec();
    // The counted API accumulates rather than replacing a caller's counter.
    let (mut old_ops, mut new_ops) = (7, 7);
    let rank = reference(&mut a, cols, &mut old_ops);
    assert_eq!(rank, rref_f2_counted(&mut b, cols, &mut new_ops));
    assert_eq!(a, b, "full RREF, including nonpivot rows and padding");
    assert!(new_ops <= old_ops);
    let mut previous = None;
    for r in 0..rank {
        let pivot = (0..cols)
            .find(|&c| b[r][c / 64] & (1 << (c % 64)) != 0)
            .unwrap();
        assert!(previous.is_none_or(|p| p < pivot));
        assert!((0..b.len()).all(|s| s == r || b[s][pivot / 64] & (1 << (pivot % 64)) == 0));
        previous = Some(pivot);
    }
    assert!(b[rank..]
        .iter()
        .all(|row| (0..cols).all(|c| row[c / 64] & (1 << (c % 64)) == 0)));
}

#[test]
fn rref_matches_reference_across_shapes_and_word_boundaries() {
    for rows in [0, 1, 2, 7, 65, 129] {
        for cols in [0usize, 1, 2, 63, 64, 65, 127, 128, 129, 257] {
            check(&vec![vec![0; cols.div_ceil(64)]; rows], cols);
            for kind in ["dense", "sparse", "deficient"] {
                for seed in [17, 937, 20260914] {
                    check(&synthetic(rows, cols, kind, seed), cols);
                }
            }
        }
    }
}

#[test]
fn rref_skips_zero_prefix_words_and_preserves_padding() {
    // Pivots start beyond a full empty word, including a skipped-column gap.
    let input = vec![vec![0, 2, u64::MAX], vec![0, 2, 6], vec![0, 0, 4]];
    check(&input, 131);
    let mut a = input.clone();
    let mut b = input;
    let (mut old, mut new) = (0, 0);
    reference(&mut a, 131, &mut old);
    rref_f2_counted(&mut b, 131, &mut new);
    assert!(new < old);
}

#[test]
fn solver_linear_tail_matches_full_rref_row_space_intersection() {
    for rows in [0, 1, 7, 65, 129] {
        for cols in [1usize, 63, 64, 65, 129, 257] {
            let low_width = cols.min(33);
            let low_start = cols - low_width;
            for kind in ["dense", "sparse", "deficient"] {
                for seed in [17, 937, 20260914] {
                    let input = synthetic(rows, cols, kind, seed);
                    let mut full = input.clone();
                    let mut candidate = input.clone();
                    let mut full_ops = 0;
                    let full_rank = rref_f2_counted(&mut full, cols, &mut full_ops);
                    let mut candidate_ops = 0;
                    let (low, low_rank) =
                        f4_solver_linear_tail(&mut candidate, cols, low_start, &mut candidate_ops);

                    let words_per_row = cols.div_ceil(64);
                    let mut flat = FlatF2Matrix {
                        data: input.iter().flatten().copied().collect(),
                        rows,
                        words: words_per_row,
                    };
                    let mut flat_ops = 0;
                    let (flat_low, flat_rank) =
                        f4_solver_linear_tail_flat(&mut flat, cols, low_start, &mut flat_ops);
                    assert_eq!(flat_rank, low_rank);
                    assert_eq!(&flat_low[..flat_rank], &low[..low_rank]);

                    let words = low_width.div_ceil(64).max(1);
                    let mut expected = Vec::new();
                    for source in full.iter().take(full_rank) {
                        let has_high = (0..low_start)
                            .any(|column| source[column / 64] & (1u64 << (column % 64)) != 0);
                        if has_high {
                            continue;
                        }
                        let mut row = vec![0u64; words];
                        for column in low_start..cols {
                            if source[column / 64] & (1u64 << (column % 64)) != 0 {
                                let target = column - low_start;
                                row[target / 64] |= 1u64 << (target % 64);
                            }
                        }
                        if row.iter().any(|&word| word != 0) {
                            expected.push(row);
                        }
                    }
                    assert_eq!(low_rank, expected.len());
                    assert_eq!(&low[..low_rank], expected.as_slice());
                }
            }
        }
    }
}

#[test]
fn masked_multiplier_schedule_matches_the_full_schedule() {
    for n_vars in 0..=12 {
        for degree in 0..=4 {
            assert_eq!(
                monomials_up_to_mask(all_variable_mask(n_vars), degree),
                monomials_up_to(n_vars, degree),
            );
        }
    }
}

#[test]
fn cached_layout_requires_exact_column_support() {
    let columns = vec![3u64, 2, 1, 0];
    let layout = F4ColumnLayout::new(columns);
    assert!(pack_rows_with_layout(&[vec![3, 0], vec![2, 1]], &layout, true).is_some());
    assert!(pack_rows_with_layout(&[vec![3, 0], vec![2]], &layout, true).is_none());
    assert!(pack_rows_with_layout(&[vec![4]], &layout, true).is_none());
}

#[test]
fn fused_flat_packing_matches_materialized_rows() {
    let n_vars = 8;
    let polynomials = vec![
        F2BoolPoly::from_monos(
            vec![
                F2BoolMono::from_mask(0b0000_0011),
                F2BoolMono::from_mask(0b0001_0100),
                F2BoolMono::var(6),
                F2BoolMono::one(),
            ],
            n_vars,
        ),
        F2BoolPoly::from_monos(
            vec![
                F2BoolMono::from_mask(0b0010_0100),
                F2BoolMono::from_mask(0b0100_1000),
                F2BoolMono::var(1),
            ],
            n_vars,
        ),
    ];
    let multiplier_mask = occurring_vars(&polynomials);
    for degree in [2, 3, 4] {
        let rows =
            macaulay_rows_monos_with_mask(&polynomials, n_vars, degree, multiplier_mask, None)
                .unwrap();
        let columns = macaulay_columns(&rows).unwrap();
        let layout = F4ColumnLayout::new(columns);
        let materialized = pack_rows_flat_with_layout(&rows, &layout, true).unwrap();
        let fused =
            pack_polynomials_flat_fused(&polynomials, n_vars, degree, multiplier_mask, &layout)
                .unwrap();
        assert_eq!(fused.rows, materialized.rows);
        assert_eq!(fused.words, materialized.words);
        assert_eq!(fused.data, materialized.data);

        let mut missing_columns = layout.columns.clone();
        missing_columns.pop();
        let missing = F4ColumnLayout::new(missing_columns);
        assert!(pack_polynomials_flat_fused(
            &polynomials,
            n_vars,
            degree,
            multiplier_mask,
            &missing,
        )
        .is_none());

        let mut extra_columns = layout.columns.clone();
        extra_columns.push(1u64 << 63);
        let extra = F4ColumnLayout::new(extra_columns);
        assert!(
            pack_polynomials_flat_fused(&polynomials, n_vars, degree, multiplier_mask, &extra,)
                .is_none()
        );
    }
}

#[test]
fn fast_column_index_matches_standard_hashing() {
    let mut columns = Vec::new();
    let mut state = 0x9e3779b97f4a7c15u64;
    for _ in 0..2048 {
        state = state.wrapping_add(0x9e3779b97f4a7c15);
        let mut value = state;
        value = (value ^ (value >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        value = (value ^ (value >> 27)).wrapping_mul(0x94d049bb133111eb);
        columns.push(value ^ (value >> 31));
    }
    columns.sort_unstable();
    columns.dedup();
    let standard: std::collections::HashMap<u64, usize> = columns
        .iter()
        .enumerate()
        .map(|(column, &monomial)| (monomial, column))
        .collect();
    let mut fast = FastColumnMap::with_capacity_and_hasher(
        columns.len(),
        std::hash::BuildHasherDefault::default(),
    );
    fast.extend(
        columns
            .iter()
            .enumerate()
            .map(|(column, &monomial)| (monomial, column)),
    );
    for monomial in &columns {
        assert_eq!(fast.get(monomial), standard.get(monomial));
    }
    assert_eq!(fast.get(&0xDEADBEEF), standard.get(&0xDEADBEEF));
}

#[test]
#[ignore = "paired release benchmark; run via research/f4_linear_algebra_20260914/run.py"]
fn paired_kernel_benchmark() {
    let contract: serde_json::Value = serde_json::from_str(include_str!(
        "../../research/f4_linear_algebra_20260914/contract.json"
    ))
    .unwrap();
    let repeats = contract["kernel"]["paired_repetitions"].as_u64().unwrap();
    for seed in [17, 937] {
        let mut cases = Vec::new();
        for spec in contract["kernel"]["synthetic"].as_array().unwrap() {
            let rows = spec[0].as_u64().unwrap() as usize;
            let cols = spec[1].as_u64().unwrap() as usize;
            let kind = spec[2].as_str().unwrap();
            let start = Instant::now();
            let matrix = synthetic(rows, cols, kind, seed);
            cases.push((
                format!("{kind}-{rows}x{cols}"),
                cols,
                matrix,
                start.elapsed().as_nanos(),
            ));
        }
        for spec in contract["kernel"]["macaulay"].as_array().unwrap() {
            let n = spec[0].as_u64().unwrap() as u32;
            let ell = spec[1].as_u64().unwrap() as usize;
            let m = spec[2].as_u64().unwrap() as usize;
            let degree = spec[3].as_u64().unwrap() as u32;
            let start = Instant::now();
            let fe = |x| F2mElement::from_biguint(&num_bigint::BigUint::from(x), n);
            let st = FieldStructure::new(
                n,
                &crate::cryptanalysis::koblitz_index_calculus::find_irreducible(n).unwrap(),
            );
            let basis: Vec<_> = (0..ell).map(|k| fe(1u64 << k)).collect();
            let sys =
                build_decomposition_system(&basis, &fe(seed % (1 << n)), &fe(1), m, &st).unwrap();
            let (columns, matrix) = build_macaulay(&sys.equations, sys.n_vars, degree).unwrap();
            cases.push((
                format!("macaulay-n{n}-l{ell}-m{m}-d{degree}"),
                columns.len(),
                matrix,
                start.elapsed().as_nanos(),
            ));
        }
        for (name, cols, input, setup_ns) in cases {
            let input_hash = blake3::hash(&serde_json::to_vec(&input).unwrap())
                .to_hex()
                .to_string();
            check(&input, cols); // correctness and one untimed warmup per kernel
            for rep in 0..repeats {
                let start = Instant::now();
                let mut a = input.clone();
                let mut b = input.clone();
                let copy_ns = start.elapsed().as_nanos();
                let mut ns = [0, 0];
                let mut ops = [0, 0];
                let mut ranks = [0, 0];
                for variant in if rep % 2 == 0 { [0, 1] } else { [1, 0] } {
                    let start = Instant::now();
                    ranks[variant] = if variant == 0 {
                        reference(black_box(&mut a), cols, &mut ops[variant])
                    } else {
                        rref_f2_counted(black_box(&mut b), cols, &mut ops[variant])
                    };
                    ns[variant] = start.elapsed().as_nanos();
                }
                let start = Instant::now();
                assert_eq!(ranks[0], ranks[1]);
                assert_eq!(a, b);
                let verification_ns = start.elapsed().as_nanos();
                println!(
                    "{}",
                    json!({"phase":"kernel","case":name,"seed":seed,"rep":rep,
                    "rows":input.len(),"cols":cols,"rank":ranks[0],"input_hash":input_hash,
                    "setup_ns":setup_ns,"copy_ns":copy_ns,"verification_ns":verification_ns,
                    "reference_ns":ns[0],"candidate_ns":ns[1],"reference_xors":ops[0],
                    "candidate_xors":ops[1],"correct":true})
                );
            }
        }
    }
}
