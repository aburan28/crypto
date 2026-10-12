// Necessary affine fibers obtained by an exact quotient of equation space.
struct EquationProjection {
    rows: [u32; 32],
    markers: [u32; 32],
    rank: usize,
}
impl EquationProjection {
    fn new(form: &SyndromeForm, k: usize) -> Self {
        assert!(k <= form.n && k <= 6);
        let mut out = Self {
            rows: [0; 32],
            markers: [0; 32],
            rank: 0,
        };
        for j in 0..k {
            for i in 0..j {
                let value = out.apply(form.quadratic[i][j]);
                if value != 0 {
                    out.rows[out.rank] = value;
                    out.markers[out.rank] = value & value.wrapping_neg();
                    out.rank += 1;
                }
            }
        }
        out
    }
    fn apply(&self, mut value: u32) -> u32 {
        for p in 0..self.rank {
            if value & self.markers[p] != 0 {
                value ^= self.rows[p];
            }
        }
        value
    }
    fn form(&self, original: &SyndromeForm) -> SyndromeForm {
        let mut out = SyndromeForm::new(original.n);
        out.constant = self.apply(original.constant);
        for i in 0..original.n {
            out.linear[i] = self.apply(original.linear[i]);
            for j in i + 1..original.n {
                let value = self.apply(original.quadratic[i][j]);
                out.quadratic[i][j] = value;
                out.quadratic[j][i] = value;
            }
        }
        out
    }
}
#[derive(Debug)]
struct AffineFiber {
    particular: Option<u32>,
    null_basis: [u32; 6],
    nullity: usize,
    rank: usize,
}
fn affine_fiber<const K: usize>(columns: &[u32; K], mut rhs: u32) -> AffineFiber {
    assert!(K <= 6);
    let mut rows = [0u32; K];
    let mut tags = [0u32; K];
    let mut markers = [0u32; K];
    let mut out = AffineFiber {
        particular: None,
        null_basis: [0; 6],
        nullity: 0,
        rank: 0,
    };
    for (j, &column) in columns.iter().enumerate() {
        let mut value = column;
        let mut tag = 1u32 << j;
        for p in 0..out.rank {
            if value & markers[p] != 0 {
                value ^= rows[p];
                tag ^= tags[p];
            }
        }
        if value == 0 {
            out.null_basis[out.nullity] = tag;
            out.nullity += 1;
        } else {
            rows[out.rank] = value;
            tags[out.rank] = tag;
            markers[out.rank] = value & value.wrapping_neg();
            out.rank += 1;
        }
    }
    let mut particular = 0;
    for p in 0..out.rank {
        if rhs & markers[p] != 0 {
            rhs ^= rows[p];
            particular ^= tags[p];
        }
    }
    if rhs == 0 {
        out.particular = Some(particular);
    }
    out
}
#[derive(Default, Clone, Debug, PartialEq, Eq)]
struct ProjectedWork {
    low_variables: usize,
    annihilated_rank: usize,
    quotient_dimension: usize,
    prefixes: u64,
    batches: u64,
    screen_rejected: u64,
    affine_queries: u64,
    affine_rejected: u64,
    rank_sum: u64,
    consistent_rank_counts: [u64; 7],
    extension_space: u64,
    extensions_checked: u64,
    original_rejected: u64,
}
impl ProjectedWork {
    fn json(&self) -> String {
        format!("{{\"low_variables\":{},\"annihilated_rank\":{},\"quotient_dimension\":{},\"prefixes\":{},\"batches\":{},\"screen_rejected\":{},\"affine_queries\":{},\"affine_rejected\":{},\"rank_sum\":{},\"consistent_rank_counts\":{:?},\"extension_space\":{},\"extensions_checked\":{},\"original_rejected\":{}}}", self.low_variables,self.annihilated_rank,self.quotient_dimension,self.prefixes,self.batches,self.screen_rejected,self.affine_queries,self.affine_rejected,self.rank_sum,self.consistent_rank_counts,self.extension_space,self.extensions_checked,self.original_rejected)
    }
}
fn projected_result(outcome: Outcome, work: ProjectedWork) -> (Solved, ProjectedWork) {
    (
        Solved {
            outcome,
            logical: Logical::default(),
            profile: Profile::default(),
            trace: 0,
        },
        work,
    )
}
fn solve_projected<const K: usize>(
    system: &System,
    n: u8,
    prefix_cap: u64,
    extension_cap: u64,
) -> (Solved, ProjectedWork) {
    let mut work = ProjectedWork {
        low_variables: K,
        ..ProjectedWork::default()
    };
    if n > 24 || K > 6 || K > n as usize {
        return projected_result(Outcome::Unknown("PROJECTED_DOMAIN"), work);
    }
    let Some(original) = SyndromeForm::from_system(system, n) else {
        return projected_result(Outcome::Unknown("PROJECTED_DOMAIN"), work);
    };
    let projection = EquationProjection::new(&original, K);
    work.annihilated_rank = projection.rank;
    work.quotient_dimension = system.len() - projection.rank;
    let projected = projection.form(&original);
    let plan = FiberPlan::new(&projected, system.len(), (1u32 << K) - 1);
    let size = 1usize << plan.outside.len().min(4);
    let high = plan.outside.len().saturating_sub(4);
    let mut cursor = GrayCursor::new(&plan.form);
    let offsets = plan.form.low_offsets();
    let mut columns: [u32; K] = std::array::from_fn(|i| plan.columns[i]);
    for step in 0..1u64 << high {
        if prefix_cap.saturating_sub(work.prefixes) < size as u64 {
            return projected_result(Outcome::Unknown("PROJECTED_PREFIX_CAP"), work);
        }
        cursor.advance(&plan.form, step);
        if step != 0 {
            let j = step.trailing_zeros() as usize;
            for (i, column) in columns.iter_mut().enumerate() {
                *column ^= plan.cross_high[j][i];
            }
        }
        let screen = if size == 16 {
            #[cfg(any(target_arch = "aarch64", target_arch = "x86_64"))]
            {
                screen_fibers_native_small::<K>(&cursor, &offsets, &columns, &plan.column_offsets)
            }
            #[cfg(not(any(target_arch = "aarch64", target_arch = "x86_64")))]
            {
                screen_fibers_scalar(&cursor, &offsets, &columns, &plan.column_offsets, size)
            }
        } else {
            screen_fibers_scalar(&cursor, &offsets, &columns, &plan.column_offsets, size)
        };
        work.prefixes += size as u64;
        work.batches += 1;
        work.screen_rejected += size as u64 - u64::from(screen.survivors.count_ones());
        let mut survivors = screen.survivors;
        while survivors != 0 {
            let y = survivors.trailing_zeros() as usize;
            survivors &= survivors - 1;
            let a: [u32; K] =
                std::array::from_fn(|i| columns[i] ^ plan.column_offsets[i][y / 4][y % 4]);
            let fiber = affine_fiber(&a, screen.values[y / 4][y % 4]);
            work.affine_queries += 1;
            work.rank_sum += fiber.rank as u64;
            let Some(mut inside) = fiber.particular else {
                work.affine_rejected += 1;
                continue;
            };
            work.consistent_rank_counts[fiber.rank] += 1;
            work.extension_space += 1u64 << fiber.nullity;
            let outside = ((step ^ (step >> 1)) << 4) | y as u64;
            for extension in 0..1u64 << fiber.nullity {
                if work.extensions_checked == extension_cap {
                    return projected_result(Outcome::Unknown("PROJECTED_EXTENSION_CAP"), work);
                }
                if extension != 0 {
                    inside ^= fiber.null_basis[extension.trailing_zeros() as usize];
                }
                let point = plan.recover(outside, inside);
                work.extensions_checked += 1;
                // Exact original equations are required even when the quotient is zero.
                if satisfies(system, point) {
                    return projected_result(Outcome::Sat(point), work);
                }
                work.original_rejected += 1;
            }
        }
    }
    projected_result(Outcome::Unsat, work)
}

#[cfg(test)]
mod projected_tests {
    use super::*;
    fn all_fiber(columns: &[u32; 3], rhs: u32) -> Vec<u32> {
        let answer = affine_fiber(columns, rhs);
        let mut out = Vec::new();
        if let Some(base) = answer.particular {
            for bits in 0..1usize << answer.nullity {
                let mut point = base;
                for j in 0..answer.nullity {
                    if bits & (1 << j) != 0 {
                        point ^= answer.null_basis[j];
                    }
                }
                out.push(point);
            }
        }
        out.sort_unstable();
        out
    }
    #[test]
    fn complete_affine_fibers_cover_every_small_column_system() {
        for code in 0..512u32 {
            let columns = [code & 7, (code >> 3) & 7, (code >> 6) & 7];
            for rhs in 0..8 {
                let direct: Vec<_> = (0..8)
                    .filter(|x| {
                        (0..3).fold(0, |v, j| if x & (1 << j) != 0 { v ^ columns[j] } else { v })
                            == rhs
                    })
                    .collect();
                assert_eq!(all_fiber(&columns, rhs), direct);
            }
        }
    }
    #[test]
    fn zero_projection_needs_nonparticular_extension_and_original_check() {
        let system = vec![vec![3, 0]];
        let (got, work) = solve_projected::<2>(&system, 2, 1, 4);
        assert_eq!(got.outcome, Outcome::Sat(3));
        assert_eq!(work.quotient_dimension, 0);
        assert_eq!(work.original_rejected, 2);
        assert_eq!(work.extensions_checked, 3);
        let (limited, _) = solve_projected::<2>(&system, 2, 1, 2);
        assert_eq!(limited.outcome, Outcome::Unknown("PROJECTED_EXTENSION_CAP"));
    }
    #[test]
    fn every_pair_of_three_variable_quadratics_matches_truth_table() {
        let monomials = [0, 1, 2, 4, 3, 5, 6];
        let poly = |code: usize| {
            monomials
                .iter()
                .enumerate()
                .filter_map(|(i, &m)| (code & (1 << i) != 0).then_some(m))
                .collect::<Vec<_>>()
        };
        for a in 0..128 {
            for b in 0..128 {
                let system = vec![poly(a), poly(b)];
                let expected = (0..8).any(|x| satisfies(&system, x));
                for got in [
                    solve_projected::<1>(&system, 3, 8, 8),
                    solve_projected::<2>(&system, 3, 8, 8),
                    solve_projected::<3>(&system, 3, 8, 8),
                ] {
                    match got.0.outcome {
                        Outcome::Sat(x) => {
                            assert!(expected && satisfies(&system, x));
                        }
                        Outcome::Unsat => {
                            assert!(!expected);
                            assert_eq!(got.1.extensions_checked, got.1.extension_space);
                        }
                        Outcome::Unknown(_) => panic!("uncapped finite control must complete"),
                    }
                }
            }
        }
    }
    #[test]
    fn projected_coefficients_match_direct_evaluation_under_all_small_restrictions() {
        let mut rng = 930172;
        for n in 2..=8 {
            for _ in 0..16 {
                let mut system = vec![Vec::new(); 7];
                for p in &mut system {
                    for m in 0..1u64 << n {
                        if m.count_ones() <= 2 && next(&mut rng) & 3 == 0 {
                            p.push(m);
                        }
                    }
                }
                let original = SyndromeForm::from_system(&system, n).unwrap();
                for k in 1..=usize::from(n).min(6) {
                    let projection = EquationProjection::new(&original, k);
                    let projected = projection.form(&original);
                    for x in 0..1u64 << n {
                        assert_eq!(projected.value(x), projection.apply(original.value(x)));
                    }
                    for j in 0..k {
                        for i in 0..j {
                            assert_eq!(projected.quadratic[i][j], 0);
                        }
                    }
                }
            }
        }
    }
    #[test]
    fn changing_rank_constants_and_caps_preserve_unknown() {
        // In y(1+z), the low coefficient disappears when z=1.
        let form = SyndromeForm::from_system(&vec![vec![1, 3]], 2).unwrap();
        let plan = FiberPlan::new(&form, 1, 1);
        assert_eq!(affine_fiber(&[plan.columns[0]], 0).rank, 1);
        assert_eq!(
            affine_fiber(&[plan.columns[0] ^ plan.column_offsets[0][0][1]], 0).rank,
            0
        );
        let inconsistent = vec![vec![0]];
        for cap in 0..8 {
            assert!(matches!(
                solve_projected::<1>(&inconsistent, 4, cap, 16).0.outcome,
                Outcome::Unknown(_)
            ));
        }
        assert_eq!(
            solve_projected::<1>(&inconsistent, 4, 8, 16).0.outcome,
            Outcome::Unsat
        );
        assert_eq!(
            solve_projected::<4>(&vec![], 4, 1, 0).0.outcome,
            Outcome::Unknown("PROJECTED_EXTENSION_CAP")
        );
        assert_eq!(
            solve_projected::<4>(&vec![vec![7]], 4, 16, 16).0.outcome,
            Outcome::Unknown("PROJECTED_DOMAIN")
        );
    }
    #[test]
    fn complete_generated_systems_match_exhaustive_search() {
        for seed in [17, 937, 7119] {
            for family in ["planted", "cross_planted", "unplanted"] {
                let (system, _) = fixture(10, seed, family);
                let exists = (0..1024).any(|x| satisfies(&system, x));
                for (got, work) in [
                    solve_projected::<4>(&system, 10, 1024, 1024),
                    solve_projected::<5>(&system, 10, 1024, 1024),
                    solve_projected::<6>(&system, 10, 1024, 1024),
                ] {
                    match got.outcome {
                        Outcome::Sat(x) => assert!(exists && satisfies(&system, x)),
                        Outcome::Unsat => assert!(!exists),
                        _ => panic!("complete expected"),
                    }
                    assert_eq!(
                        work.affine_queries - work.affine_rejected,
                        work.consistent_rank_counts.iter().sum()
                    );
                }
            }
        }
    }
}
