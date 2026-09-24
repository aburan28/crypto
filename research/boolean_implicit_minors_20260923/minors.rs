// Maximal augmented minors are necessary consistency conditions, including at
// rank drops. Boolean products use monomial union and coefficient parity.
const MINOR_K: usize = 4;
const MINOR_COUNT: usize = 4;
type MinorMatrix = Vec<Vec<Vec<u16>>>;
#[derive(Default, Clone, Debug, PartialEq, Eq)]
struct MinorWork {
    annihilated_rank: usize,
    quotient_dimension: usize,
    product_toggles: u64,
    cancelled_toggles: u64,
    support_sum: u64,
    dp_states: u64,
    max_support: usize,
    dp_payload_peak_bytes: usize,
    truth_word_xors: u64,
    numeric_determinants: u64,
    prefixes: u64,
    affine_queries: u64,
    affine_rejected: u64,
    screen_rejected: u64,
    rank_sum: u64,
    degrees: Vec<Option<u32>>,
    supports: Vec<usize>,
}
impl MinorWork {
    fn json(&self) -> String {
        let degrees = self
            .degrees
            .iter()
            .map(|x| x.map_or("null".into(), |v| v.to_string()))
            .collect::<Vec<_>>()
            .join(",");
        format!("{{\"annihilated_rank\":{},\"quotient_dimension\":{},\"product_toggles\":{},\"cancelled_toggles\":{},\"support_sum\":{},\"dp_states\":{},\"max_support\":{},\"dp_payload_peak_bytes\":{},\"truth_word_xors\":{},\"numeric_determinants\":{},\"prefixes\":{},\"affine_queries\":{},\"affine_rejected\":{},\"screen_rejected\":{},\"rank_sum\":{},\"degrees\":[{}],\"supports\":{:?}}}",self.annihilated_rank,self.quotient_dimension,self.product_toggles,self.cancelled_toggles,self.support_sum,self.dp_states,self.max_support,self.dp_payload_peak_bytes,self.truth_word_xors,self.numeric_determinants,self.prefixes,self.affine_queries,self.affine_rejected,self.screen_rejected,self.rank_sum,degrees,self.supports)
    }
}
#[derive(Clone, Debug, PartialEq, Eq)]
struct MinorOutput {
    accepted: Vec<u64>,
    minors: Vec<Vec<u64>>,
}
struct MinorRun {
    output: Option<MinorOutput>,
    reason: Option<&'static str>,
    work: MinorWork,
    setup_ns: u128,
    construction_ns: u128,
    evaluation_ns: u128,
}
fn minor_words(h: usize) -> usize {
    (1usize << h).div_ceil(64)
}
fn minor_full_mask(h: usize) -> Vec<u64> {
    let size = 1usize << h;
    let mut mask = vec![u64::MAX; minor_words(h)];
    if size < 64 {
        mask[0] = (1u64 << size) - 1;
    }
    mask
}
fn minor_accept(truth: &[Vec<u64>], h: usize) -> Vec<u64> {
    let mut mask = minor_full_mask(h);
    for t in truth {
        for (dst, &word) in mask.iter_mut().zip(t) {
            *dst &= !word;
        }
    }
    mask
}
fn minor_rows(projection: &EquationProjection, equations: usize) -> Option<Vec<[usize; 5]>> {
    let pivots = projection.markers[..projection.rank]
        .iter()
        .fold(0, |v, &p| v | p);
    let available: Vec<_> = (0..equations)
        .filter(|&e| pivots & (1u32 << e) == 0)
        .collect();
    if available.len() < 5 {
        return None;
    }
    Some(
        (0..MINOR_COUNT)
            .map(|offset| std::array::from_fn(|j| available[(offset + j) % available.len()]))
            .collect(),
    )
}
fn minor_matrix(form: &SyndromeForm, rows: &[usize; 5]) -> MinorMatrix {
    let h = form.n - MINOR_K;
    rows.iter()
        .map(|&e| {
            let mut row = Vec::new();
            for i in 0..MINOR_K {
                let mut poly = Vec::new();
                if form.linear[i] & (1u32 << e) != 0 {
                    poly.push(0);
                }
                for j in 0..h {
                    if form.quadratic[i][j + MINOR_K] & (1u32 << e) != 0 {
                        poly.push(1u16 << j);
                    }
                }
                row.push(poly);
            }
            let mut rhs = Vec::new();
            if form.constant & (1u32 << e) != 0 {
                rhs.push(0);
            }
            for j in 0..h {
                if form.linear[j + MINOR_K] & (1u32 << e) != 0 {
                    rhs.push(1u16 << j);
                }
                for i in 0..j {
                    if form.quadratic[i + MINOR_K][j + MINOR_K] & (1u32 << e) != 0 {
                        rhs.push((1u16 << i) | (1u16 << j));
                    }
                }
            }
            rhs.sort_unstable();
            row.push(rhs);
            row
        })
        .collect()
}
fn minor_compile(
    matrix: &MinorMatrix,
    h: usize,
    cap: u64,
    work: &mut MinorWork,
) -> Result<Vec<u64>, &'static str> {
    assert!(h <= 12 && matrix.len() == 5 && matrix.iter().all(|r| r.len() == 5));
    let words = minor_words(h);
    let mut dp = vec![vec![0u64; words]; 32];
    work.dp_payload_peak_bytes = work.dp_payload_peak_bytes.max(32 * words * 8);
    dp[0][0] = 1;
    for mask in 1usize..32 {
        let column = mask.count_ones() as usize - 1;
        let mut support = 0usize;
        for (row, entries) in matrix.iter().enumerate() {
            if mask & (1 << row) == 0 {
                continue;
            }
            let source = mask ^ (1 << row);
            for word in 0..words {
                let mut bits = dp[source][word];
                while bits != 0 {
                    let monomial = word * 64 + bits.trailing_zeros() as usize;
                    bits &= bits - 1;
                    for &term in &entries[column] {
                        if work.product_toggles == cap {
                            return Err("PRODUCT_CAP");
                        }
                        assert!(usize::from(term) < 1usize << h);
                        let target = monomial | usize::from(term);
                        let bit = 1u64 << (target % 64);
                        let cancelled = dp[mask][target / 64] & bit != 0;
                        dp[mask][target / 64] ^= bit;
                        work.product_toggles += 1;
                        if cancelled {
                            work.cancelled_toggles += 1;
                            work.support_sum -= 1;
                            support -= 1;
                        } else {
                            work.support_sum += 1;
                            support += 1;
                        }
                        work.max_support = work.max_support.max(support);
                    }
                }
            }
        }
        work.dp_states += 1;
    }
    Ok(std::mem::take(&mut dp[31]))
}
fn minor_truth(mut coefficients: Vec<u64>, h: usize, work: &mut MinorWork) -> Vec<u64> {
    // Boolean zeta/Mobius transform, packed within each word for the first six
    // variables and between words thereafter. It is its own inverse over F2.
    for j in 0..h {
        if j < 6 {
            let shift = 1usize << j;
            let mut upper = 0u64;
            for p in 0..64 {
                if p & shift != 0 {
                    upper |= 1u64 << p;
                }
            }
            for word in &mut coefficients {
                *word ^= (*word << shift) & upper;
                work.truth_word_xors += 1;
            }
        } else {
            let stride = 1usize << (j - 6);
            for start in (0..coefficients.len()).step_by(2 * stride) {
                for i in 0..stride {
                    coefficients[start + stride + i] ^= coefficients[start + i];
                    work.truth_word_xors += 1;
                }
            }
        }
    }
    if h < 6 {
        coefficients[0] &= (1u64 << (1usize << h)) - 1;
    }
    coefficients
}
fn minor_degree(coefficients: &[u64]) -> Option<u32> {
    let mut degree = None;
    for (w, &word) in coefficients.iter().enumerate() {
        let mut bits = word;
        while bits != 0 {
            let m = w * 64 + bits.trailing_zeros() as usize;
            bits &= bits - 1;
            degree = Some(degree.unwrap_or(0).max(m.count_ones()));
        }
    }
    degree
}
fn minor_det5(mut rows: [u8; 5]) -> bool {
    for column in 0..5 {
        let Some(pivot) = (column..5).find(|&r| rows[r] & (1 << column) != 0) else {
            return false;
        };
        rows.swap(column, pivot);
        for r in column + 1..5 {
            if rows[r] & (1 << column) != 0 {
                rows[r] ^= rows[column];
            }
        }
    }
    true
}
fn minor_numeric_rows(columns: &[u32; 4], rhs: u32, rows: &[usize; 5]) -> [u8; 5] {
    std::array::from_fn(|r| {
        let mut value = (((rhs >> rows[r]) & 1) as u8) << 4;
        for (i, &column) in columns.iter().enumerate() {
            value |= (((column >> rows[r]) & 1) as u8) << i;
        }
        value
    })
}
fn minor_scan(
    form: &SyndromeForm,
    equations: usize,
    rows: &[[usize; 5]],
    affine: bool,
    work: &mut MinorWork,
) -> MinorOutput {
    let plan = FiberPlan::new(form, equations, (1u32 << MINOR_K) - 1);
    let h = plan.outside.len();
    let size = 1usize << h.min(4);
    let mut cursor = GrayCursor::new(&plan.form);
    let offsets = plan.form.low_offsets();
    let mut columns: [u32; 4] = std::array::from_fn(|i| plan.columns[i]);
    let mut out = MinorOutput {
        accepted: vec![0; minor_words(h)],
        minors: if affine {
            vec![]
        } else {
            vec![vec![0; minor_words(h)]; rows.len()]
        },
    };
    for step in 0..1u64 << h.saturating_sub(4) {
        cursor.advance(&plan.form, step);
        if step != 0 {
            let j = step.trailing_zeros() as usize;
            for (i, c) in columns.iter_mut().enumerate() {
                *c ^= plan.cross_high[j][i];
            }
        }
        let screen = if affine {
            Some(if size == 16 {
                screen_fibers_native(&cursor, &offsets, &columns, &plan.column_offsets)
            } else {
                screen_fibers_scalar(&cursor, &offsets, &columns, &plan.column_offsets, size)
            })
        } else {
            None
        };
        work.prefixes += size as u64;
        for y in 0..size {
            if screen.as_ref().is_some_and(|s| s.survivors & (1 << y) == 0) {
                work.screen_rejected += 1;
                continue;
            }
            let a: [u32; 4] =
                std::array::from_fn(|i| columns[i] ^ plan.column_offsets[i][y / 4][y % 4]);
            let b = screen.as_ref().map_or_else(
                || {
                    let mut value = cursor.constant ^ offsets[y / 4][y % 4];
                    for j in 0..4 {
                        if y & (1 << j) != 0 {
                            value ^= cursor.linear[j];
                        }
                    }
                    value
                },
                |s| s.values[y / 4][y % 4],
            );
            let point = (((step ^ (step >> 1)) << 4) | y as u64) as usize;
            if affine {
                let solved = fiber_columns(&a, b);
                work.affine_queries += 1;
                work.rank_sum += u64::from(solved.rank);
                if solved.model.is_some() {
                    out.accepted[point / 64] |= 1u64 << (point % 64);
                } else {
                    work.affine_rejected += 1;
                }
            } else {
                for (i, r) in rows.iter().enumerate() {
                    work.numeric_determinants += 1;
                    if minor_det5(minor_numeric_rows(&a, b, r)) {
                        out.minors[i][point / 64] |= 1u64 << (point % 64);
                    }
                }
            }
        }
    }
    if !affine {
        out.accepted = minor_accept(&out.minors, h);
    }
    out
}
fn minor_build(system: &System, n: u8, arm: &str, cap: u64) -> MinorRun {
    assert!((MINOR_K..=16).contains(&(n as usize)));
    let mut run = MinorRun {
        output: None,
        reason: None,
        work: MinorWork::default(),
        setup_ns: 0,
        construction_ns: 0,
        evaluation_ns: 0,
    };
    let tick = Instant::now();
    let original = SyndromeForm::from_system(system, n).unwrap();
    let projection = EquationProjection::new(&original, MINOR_K);
    run.work.annihilated_rank = projection.rank;
    run.work.quotient_dimension = system.len() - projection.rank;
    let Some(rows) = minor_rows(&projection, system.len()) else {
        run.reason = Some("ROW_DOMAIN");
        return run;
    };
    let form = projection.form(&original);
    let h = usize::from(n) - MINOR_K;
    run.setup_ns = tick.elapsed().as_nanos();
    let tick = Instant::now();
    if arm == "symbolic_minor4" {
        let mut coefficients = Vec::new();
        for row in &rows {
            let matrix = minor_matrix(&form, row);
            let poly = match minor_compile(&matrix, h, cap, &mut run.work) {
                Ok(p) => p,
                Err(reason) => {
                    run.reason = Some(reason);
                    run.construction_ns = tick.elapsed().as_nanos();
                    return run;
                }
            };
            run.work
                .supports
                .push(poly.iter().map(|w| w.count_ones() as usize).sum());
            run.work.degrees.push(minor_degree(&poly));
            coefficients.push(poly);
        }
        run.construction_ns = tick.elapsed().as_nanos();
        let tick = Instant::now();
        let minors: Vec<_> = coefficients
            .into_iter()
            .map(|p| minor_truth(p, h, &mut run.work))
            .collect();
        run.output = Some(MinorOutput {
            accepted: minor_accept(&minors, h),
            minors,
        });
        run.evaluation_ns = tick.elapsed().as_nanos();
    } else {
        assert!(arm == "numeric_minor4" || arm == "affine_filter");
        run.output = Some(minor_scan(
            &form,
            system.len(),
            &rows,
            arm == "affine_filter",
            &mut run.work,
        ));
        run.evaluation_ns = tick.elapsed().as_nanos();
    }
    run
}
struct MinorOracle {
    minor: MinorOutput,
    affine: MinorOutput,
    rows: Vec<[usize; 5]>,
    original_solutions: u64,
}
fn minor_oracle(system: &System, n: u8) -> MinorOracle {
    let original = SyndromeForm::from_system(system, n).unwrap();
    let projection = EquationProjection::new(&original, MINOR_K);
    let rows = minor_rows(&projection, system.len()).unwrap();
    let h = usize::from(n) - MINOR_K;
    let mut minors = vec![vec![0; minor_words(h)]; rows.len()];
    let mut affine = MinorOutput {
        accepted: vec![0; minor_words(h)],
        minors: vec![],
    };
    for z in 0..1usize << h {
        let point = (z as u64) << MINOR_K;
        let rhs = projection.apply(original.value(point));
        let columns: [u32; 4] =
            std::array::from_fn(|i| projection.apply(original.value(point | (1u64 << i))) ^ rhs);
        for (r, selected) in rows.iter().enumerate() {
            if minor_det5(minor_numeric_rows(&columns, rhs, selected)) {
                minors[r][z / 64] |= 1u64 << (z % 64);
            }
        }
        if fiber_rows(&columns, rhs, system.len()).model.is_some() {
            affine.accepted[z / 64] |= 1u64 << (z % 64);
        }
    }
    let minor = MinorOutput {
        accepted: minor_accept(&minors, h),
        minors,
    };
    for (&a, &b) in affine.accepted.iter().zip(&minor.accepted) {
        assert_eq!(a & !b, 0);
    }
    let mut original_solutions = 0;
    for point in 0..1u64 << n {
        if satisfies(system, point) {
            original_solutions += 1;
            let outside = (point >> MINOR_K) as usize;
            let bit = 1u64 << (outside % 64);
            assert_ne!(minor.accepted[outside / 64] & bit, 0);
            assert_ne!(affine.accepted[outside / 64] & bit, 0);
        }
    }
    MinorOracle {
        minor,
        affine,
        rows,
        original_solutions,
    }
}

#[cfg(test)]
mod minor_tests {
    use super::*;
    #[test]
    fn numeric_elimination_matches_explicit_permutation_parity() {
        fn parity(rows: &[u8; 5], row: usize, used: u8) -> bool {
            if row == 5 {
                return true;
            }
            let mut result = false;
            for column in 0..5 {
                let bit = 1 << column;
                if used & bit == 0 && rows[row] & bit != 0 {
                    result ^= parity(rows, row + 1, used | bit);
                }
            }
            result
        }
        let mut seed = 77109;
        for _ in 0..4096 {
            let rows = std::array::from_fn(|_| (next(&mut seed) & 31) as u8);
            assert_eq!(minor_det5(rows), parity(&rows, 0, 0));
        }
    }
    fn eval(poly: &[u16], point: usize) -> bool {
        poly.iter()
            .filter(|&&m| point & usize::from(m) == usize::from(m))
            .count()
            % 2
            == 1
    }
    #[test]
    fn packed_truth_transform_matches_every_small_polynomial_and_is_involutive() {
        for h in 0..=3 {
            for code in 0u64..1u64 << (1usize << h) {
                let truth = minor_truth(vec![code], h, &mut MinorWork::default());
                for x in 0..1usize << h {
                    let mut value = 0;
                    for m in 0..1usize << h {
                        if x & m == m {
                            value ^= (code >> m) & 1;
                        }
                    }
                    assert_eq!((truth[x / 64] >> (x % 64)) & 1, value);
                }
                assert_eq!(minor_truth(truth, h, &mut MinorWork::default()), vec![code]);
            }
        }
        let mut seed = 8419;
        for h in 6..=12 {
            let coefficients: Vec<_> = (0..minor_words(h)).map(|_| next(&mut seed)).collect();
            let truth = minor_truth(coefficients.clone(), h, &mut MinorWork::default());
            for x in 0..1usize << h {
                let mut value = 0;
                let mut sub = x;
                loop {
                    value ^= (coefficients[sub / 64] >> (sub % 64)) & 1;
                    if sub == 0 {
                        break;
                    }
                    sub = (sub - 1) & x;
                }
                assert_eq!((truth[x / 64] >> (x % 64)) & 1, value);
            }
            assert_eq!(
                minor_truth(truth, h, &mut MinorWork::default()),
                coefficients
            );
        }
    }
    #[test]
    fn symbolic_determinants_match_numeric_elimination_with_cancellations() {
        let mut seed = 311907;
        for h in 0..=5 {
            for _ in 0..32 {
                let matrix: MinorMatrix = (0..5)
                    .map(|_| {
                        (0..5)
                            .map(|_| {
                                (0..1u16 << h)
                                    .filter(|&m| m.count_ones() <= 2 && next(&mut seed) & 3 == 0)
                                    .collect()
                            })
                            .collect()
                    })
                    .collect();
                let mut work = MinorWork::default();
                let polynomial = minor_compile(&matrix, h, 50000000, &mut work).unwrap();
                let truth = minor_truth(polynomial, h, &mut work);
                for x in 0..1usize << h {
                    let numeric = std::array::from_fn(|r| {
                        (0..5).fold(0, |v, c| v | ((eval(&matrix[r][c], x) as u8) << c))
                    });
                    assert_eq!(
                        (truth[x / 64] >> (x % 64)) & 1,
                        u64::from(minor_det5(numeric))
                    );
                }
                assert_eq!(
                    work.product_toggles - 2 * work.cancelled_toggles,
                    work.support_sum
                );
            }
        }
    }
    #[test]
    fn duplicate_rows_degree_drops_and_constant_minors_remain_explicit() {
        let mut identity = vec![vec![vec![]; 5]; 5];
        for (r, row) in identity.iter_mut().enumerate() {
            row[r] = vec![0];
        }
        let one = minor_compile(&identity, 3, 1000, &mut MinorWork::default()).unwrap();
        assert_eq!(one, vec![1]);
        assert_eq!(minor_degree(&one), Some(0));
        identity[4] = identity[3].clone();
        let zero = minor_compile(&identity, 3, 1000, &mut MinorWork::default()).unwrap();
        assert_eq!(zero, vec![0]);
        assert_eq!(minor_degree(&zero), None);
        let mut diagonal = vec![vec![vec![]; 5]; 5];
        for (r, row) in diagonal.iter_mut().enumerate() {
            row[r] = vec![1];
        }
        let reduced = minor_compile(&diagonal, 3, 1000, &mut MinorWork::default()).unwrap();
        assert_eq!(reduced, vec![2]);
        assert_eq!(minor_degree(&reduced), Some(1));
    }
    #[test]
    fn rank_deficiency_makes_zero_minors_insufficient() {
        let columns = [0u32; 4];
        let rhs = 1;
        assert!(!minor_det5(minor_numeric_rows(
            &columns,
            rhs,
            &[0, 1, 2, 3, 4]
        )));
        assert!(fiber_rows(&columns, rhs, 5).model.is_none());
    }
    #[test]
    fn capped_compilation_is_not_a_small_or_zero_polynomial() {
        let matrix = vec![vec![vec![0, 1, 2, 3]; 5]; 5];
        let mut work = MinorWork::default();
        assert_eq!(minor_compile(&matrix, 2, 3, &mut work), Err("PRODUCT_CAP"));
        assert_eq!(work.product_toggles, 3);
        assert_eq!(
            work.product_toggles - 2 * work.cancelled_toggles,
            work.support_sum
        );
        assert!(work.dp_states < 31);
    }
    #[test]
    fn fixed_row_policy_and_all_filters_match_direct_original_equations() {
        for seed in [17, 937] {
            for family in ["planted", "cross_planted", "unplanted"] {
                let (system, _) = fixture(12, seed, family);
                let oracle = minor_oracle(&system, 12);
                for arm in ["numeric_minor4", "symbolic_minor4", "affine_filter"] {
                    let got = minor_build(&system, 12, arm, 50000000);
                    assert!(got.reason.is_none());
                    assert_eq!(
                        got.output.as_ref().unwrap(),
                        if arm == "affine_filter" {
                            &oracle.affine
                        } else {
                            &oracle.minor
                        }
                    );
                    if arm == "symbolic_minor4" {
                        assert!(got.work.degrees.iter().all(|d| d.is_none_or(|v| v <= 6)));
                    }
                }
                for row in oracle.rows {
                    assert_eq!(
                        row.iter().collect::<std::collections::BTreeSet<_>>().len(),
                        5
                    );
                }
            }
        }
    }
}
