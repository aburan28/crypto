// Exact affine fibers in an independent coordinate set of the union graph.
fn fiber_graph(form: &SyndromeForm) -> [u32; MAX] {
    let mut graph = [0; MAX];
    for j in 0..form.n {
        for i in 0..j {
            if form.quadratic[i][j] != 0 {
                graph[i] |= 1 << j;
                graph[j] |= 1 << i;
            }
        }
    }
    graph
}
fn better_independent(a: u32, b: u32) -> u32 {
    if a.count_ones() > b.count_ones() || (a.count_ones() == b.count_ones() && a < b) {
        a
    } else {
        b
    }
}
fn maximum_independent(graph: &[u32; MAX], n: usize) -> (u32, u64) {
    assert!(n <= 24);
    if graph[..n].iter().all(|&x| x == 0) {
        return ((1u32 << n) - 1, 1);
    }
    let left = n / 2;
    let right = n - left;
    let mut best = vec![0u32; 1 << right];
    for mask in 1usize..1 << right {
        let bit = mask & mask.wrapping_neg();
        let j = bit.trailing_zeros() as usize;
        let rest = mask ^ bit;
        best[mask] = better_independent(
            best[rest],
            bit as u32 | best[rest & !(graph[left + j] >> left) as usize],
        );
    }
    let mut answer = 0;
    'subset: for mask in 0u32..1 << left {
        let mut bits = mask;
        let mut allowed = (1u32 << right) - 1;
        while bits != 0 {
            let i = bits.trailing_zeros() as usize;
            bits &= bits - 1;
            if graph[i] & mask != 0 {
                continue 'subset;
            }
            allowed &= !(graph[i] >> left);
        }
        answer = better_independent(answer, mask | (best[allowed as usize] << left));
    }
    (answer, (1 << left) + (1 << right))
}
struct FiberPlan {
    n: usize,
    equations: usize,
    selected: u32,
    low: Vec<usize>,
    outside: Vec<usize>,
    form: SyndromeForm,
    columns: Vec<u32>,
    column_offsets: Vec<[[u32; 4]; 4]>,
    cross_high: Vec<Vec<u32>>,
}
impl FiberPlan {
    fn new(original: &SyndromeForm, equations: usize, selected: u32) -> Self {
        assert!(selected < 1u32 << original.n && equations <= 32);
        let low: Vec<_> = (0..original.n)
            .filter(|&j| selected & (1 << j) != 0)
            .collect();
        let outside: Vec<_> = (0..original.n)
            .filter(|&j| selected & (1 << j) == 0)
            .collect();
        for &i in &low {
            for &j in &low {
                assert_eq!(original.quadratic[i][j], 0);
            }
        }
        let mut form = SyndromeForm::new(outside.len());
        form.constant = original.constant;
        for (i, &a) in outside.iter().enumerate() {
            form.linear[i] = original.linear[a];
            for (j, &b) in outside.iter().enumerate() {
                form.quadratic[i][j] = original.quadratic[a][b];
            }
        }
        let mut column_offsets = vec![[[0; 4]; 4]; low.len()];
        for (i, &a) in low.iter().enumerate() {
            for y in 0..16 {
                for j in 0..outside.len().min(4) {
                    if y & (1 << j) != 0 {
                        column_offsets[i][y / 4][y % 4] ^= original.quadratic[a][outside[j]];
                    }
                }
            }
        }
        let columns = low.iter().map(|&j| original.linear[j]).collect();
        let cross_high = outside
            .iter()
            .skip(4)
            .map(|&j| low.iter().map(|&i| original.quadratic[i][j]).collect())
            .collect();
        Self {
            n: original.n,
            equations,
            selected,
            low,
            outside,
            form,
            columns,
            column_offsets,
            cross_high,
        }
    }
    fn recover(&self, outside_point: u64, inside: u32) -> u64 {
        let mut point = 0;
        for (j, &label) in self.low.iter().enumerate() {
            point |= u64::from((inside >> j) & 1) << label;
        }
        for (j, &label) in self.outside.iter().enumerate() {
            point |= ((outside_point >> j) & 1) << label;
        }
        point
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct LinearAnswer {
    model: Option<u32>,
    rank: u32,
}
fn fiber_rows(columns: &[u32], rhs: u32, equations: usize) -> LinearAnswer {
    let k = columns.len();
    let constant = 1u32 << k;
    let mut basis = [0u32; MAX];
    let mut rank = 0;
    let mut inconsistent = false;
    for e in 0..equations {
        let mut row = ((rhs >> e) & 1) << k;
        for (j, &a) in columns.iter().enumerate() {
            row |= ((a >> e) & 1) << j;
        }
        loop {
            let variables = row & (constant - 1);
            if variables == 0 {
                inconsistent |= row != 0;
                break;
            }
            let p = variables.trailing_zeros() as usize;
            if basis[p] == 0 {
                basis[p] = row;
                rank += 1;
                break;
            }
            row ^= basis[p];
        }
    }
    if inconsistent {
        return LinearAnswer { model: None, rank };
    }
    let mut model = 0u32;
    for p in (0..k).rev() {
        let row = basis[p];
        if row != 0 && ((row & model).count_ones() & 1) ^ ((row >> k) & 1) != 0 {
            model |= 1 << p;
        }
    }
    LinearAnswer {
        model: Some(model),
        rank,
    }
}
fn fiber_columns_small<const K: usize>(columns: &[u32], mut rhs: u32) -> LinearAnswer {
    // Greedy independent columns in low-to-high variable order give the
    // smallest binary witness: every omitted column uses earlier columns.
    let mut basis = [0u32; K];
    let mut tags = [0u32; K];
    let mut markers = [0u32; K];
    let mut rank = 0;
    for (j, &column) in columns.iter().enumerate() {
        let mut value = column;
        let mut tag = 1u32 << j;
        for p in 0..rank {
            if value & markers[p] != 0 {
                value ^= basis[p];
                tag ^= tags[p];
            }
        }
        if value != 0 {
            basis[rank] = value;
            tags[rank] = tag;
            markers[rank] = value & value.wrapping_neg();
            rank += 1;
        }
    }
    let mut model = 0;
    for p in 0..rank {
        if rhs & markers[p] != 0 {
            rhs ^= basis[p];
            model ^= tags[p];
        }
    }
    LinearAnswer {
        model: if rhs == 0 { Some(model) } else { None },
        rank: rank as u32,
    }
}
fn fiber_columns(columns: &[u32], rhs: u32) -> LinearAnswer {
    match columns.len() {
        0 => fiber_columns_small::<0>(columns, rhs),
        1 => fiber_columns_small::<1>(columns, rhs),
        2 => fiber_columns_small::<2>(columns, rhs),
        3 => fiber_columns_small::<3>(columns, rhs),
        4 => fiber_columns_small::<4>(columns, rhs),
        5 => fiber_columns_small::<5>(columns, rhs),
        _ => fiber_columns_small::<MAX>(columns, rhs),
    }
}
const FIBER_SHIFTS: [u32; 7] = [1, 2, 3, 4, 5, 7, 8];
struct FiberScreen {
    survivors: u16,
    values: [[u32; 4]; 4],
    rounds: u64,
}
fn screen_fibers_scalar(
    cursor: &GrayCursor,
    offsets: &[[u32; 4]; 4],
    columns: &[u32],
    images: &[[[u32; 4]; 4]],
    size: usize,
) -> FiberScreen {
    let mut out = FiberScreen {
        survivors: 0,
        values: [[0; 4]; 4],
        rounds: 0,
    };
    for g in 0..size.div_ceil(4) {
        let width = (size - g * 4).min(4);
        let mut values = [[0u32; MAX]; 4];
        let mut bad = [0u32; 4];
        let mut mask = 0u16;
        for y in 0..width {
            let low = g * 4 + y;
            let mut b = cursor.constant ^ offsets[g][y];
            for j in 0..4 {
                if low & (1 << j) != 0 {
                    b ^= cursor.linear[j];
                }
            }
            let mut union = 0;
            for (i, &a) in columns.iter().enumerate() {
                values[y][i] = a ^ images[i][g][y];
                union |= values[y][i];
            }
            out.values[g][y] = b;
            bad[y] = b & !union;
            if bad[y] == 0 {
                mask |= 1 << y;
            }
        }
        out.rounds += width as u64;
        for shift in FIBER_SHIFTS {
            if mask == 0 {
                break;
            }
            mask = 0;
            for y in 0..width {
                let b = out.values[g][y];
                let mut union = 0;
                for &a in &values[y][..columns.len()] {
                    union |= a ^ (a >> shift);
                }
                bad[y] |= (b ^ (b >> shift)) & !union;
                if bad[y] == 0 {
                    mask |= 1 << y;
                }
            }
            out.rounds += width as u64;
        }
        out.survivors |= mask << (g * 4);
    }
    out
}
#[cfg(target_arch = "aarch64")]
fn screen_fibers_native_small<const K: usize>(
    cursor: &GrayCursor,
    offsets: &[[u32; 4]; 4],
    columns: &[u32],
    images: &[[[u32; 4]; 4]],
) -> FiberScreen {
    unsafe {
        use std::arch::aarch64::*;
        let c = cursor.constant;
        let l = &cursor.linear;
        let base = [c, c ^ l[0], c ^ l[1], c ^ l[0] ^ l[1]];
        let base = vld1q_u32(base.as_ptr());
        let weights = [1u32, 2, 4, 8];
        let weights = vld1q_u32(weights.as_ptr());
        let zero = vdupq_n_u32(0);
        let mut out = FiberScreen {
            survivors: 0,
            values: [[0; 4]; 4],
            rounds: 0,
        };
        for g in 0..4 {
            let delta = (if g & 1 != 0 { l[2] } else { 0 }) ^ (if g & 2 != 0 { l[3] } else { 0 });
            let b = veorq_u32(
                veorq_u32(base, vdupq_n_u32(delta)),
                vld1q_u32(offsets[g].as_ptr()),
            );
            let a: [uint32x4_t; K] = std::array::from_fn(|i| {
                veorq_u32(vdupq_n_u32(columns[i]), vld1q_u32(images[i][g].as_ptr()))
            });
            let mut union = zero;
            for &v in &a {
                union = vorrq_u32(union, v);
            }
            let mut bad = vbicq_u32(b, union);
            let mut mask = vaddvq_u32(vandq_u32(vceqq_u32(bad, zero), weights));
            out.rounds += 4;
            macro_rules! pair {
                ($s:literal) => {
                    if mask != 0 {
                        let mut union = zero;
                        for &v in &a {
                            union = vorrq_u32(union, veorq_u32(v, vshrq_n_u32::<$s>(v)));
                        }
                        bad = vorrq_u32(bad, vbicq_u32(veorq_u32(b, vshrq_n_u32::<$s>(b)), union));
                        mask = vaddvq_u32(vandq_u32(vceqq_u32(bad, zero), weights));
                        out.rounds += 4;
                    }
                };
            }
            pair!(1);
            pair!(2);
            pair!(3);
            pair!(4);
            pair!(5);
            pair!(7);
            pair!(8);
            out.survivors |= (mask as u16) << (g * 4);
            vst1q_u32(out.values[g].as_mut_ptr(), b);
        }
        out
    }
}
#[cfg(target_arch = "x86_64")]
fn screen_fibers_native_small<const K: usize>(
    cursor: &GrayCursor,
    offsets: &[[u32; 4]; 4],
    columns: &[u32],
    images: &[[[u32; 4]; 4]],
) -> FiberScreen {
    unsafe {
        use std::arch::x86_64::*;
        let c = cursor.constant;
        let l = &cursor.linear;
        let base = [c, c ^ l[0], c ^ l[1], c ^ l[0] ^ l[1]];
        let base = _mm_loadu_si128(base.as_ptr().cast());
        let zero = _mm_setzero_si128();
        let mut out = FiberScreen {
            survivors: 0,
            values: [[0; 4]; 4],
            rounds: 0,
        };
        for g in 0..4 {
            let delta = (if g & 1 != 0 { l[2] } else { 0 }) ^ (if g & 2 != 0 { l[3] } else { 0 });
            let b = _mm_xor_si128(
                _mm_xor_si128(base, _mm_set1_epi32(delta as i32)),
                _mm_loadu_si128(offsets[g].as_ptr().cast()),
            );
            let a: [__m128i; K] = std::array::from_fn(|i| {
                _mm_xor_si128(
                    _mm_set1_epi32(columns[i] as i32),
                    _mm_loadu_si128(images[i][g].as_ptr().cast()),
                )
            });
            let mut union = zero;
            for &v in &a {
                union = _mm_or_si128(union, v);
            }
            let mut bad = _mm_andnot_si128(union, b);
            let mut mask = _mm_movemask_ps(_mm_castsi128_ps(_mm_cmpeq_epi32(bad, zero)));
            out.rounds += 4;
            macro_rules! pair {
                ($s:literal) => {
                    if mask != 0 {
                        let mut union = zero;
                        for &v in &a {
                            union = _mm_or_si128(union, _mm_xor_si128(v, _mm_srli_epi32::<$s>(v)));
                        }
                        bad = _mm_or_si128(
                            bad,
                            _mm_andnot_si128(union, _mm_xor_si128(b, _mm_srli_epi32::<$s>(b))),
                        );
                        mask = _mm_movemask_ps(_mm_castsi128_ps(_mm_cmpeq_epi32(bad, zero)));
                        out.rounds += 4;
                    }
                };
            }
            pair!(1);
            pair!(2);
            pair!(3);
            pair!(4);
            pair!(5);
            pair!(7);
            pair!(8);
            out.survivors |= (mask as u16) << (g * 4);
            _mm_storeu_si128(out.values[g].as_mut_ptr().cast(), b);
        }
        out
    }
}
#[cfg(any(target_arch = "aarch64", target_arch = "x86_64"))]
fn screen_fibers_native(
    cursor: &GrayCursor,
    offsets: &[[u32; 4]; 4],
    columns: &[u32],
    images: &[[[u32; 4]; 4]],
) -> FiberScreen {
    match columns.len() {
        0 => screen_fibers_native_small::<0>(cursor, offsets, columns, images),
        1 => screen_fibers_native_small::<1>(cursor, offsets, columns, images),
        2 => screen_fibers_native_small::<2>(cursor, offsets, columns, images),
        3 => screen_fibers_native_small::<3>(cursor, offsets, columns, images),
        4 => screen_fibers_native_small::<4>(cursor, offsets, columns, images),
        5 => screen_fibers_native_small::<5>(cursor, offsets, columns, images),
        _ => screen_fibers_scalar(cursor, offsets, columns, images, 16),
    }
}
#[cfg(not(any(target_arch = "aarch64", target_arch = "x86_64")))]
fn screen_fibers_native(
    cursor: &GrayCursor,
    offsets: &[[u32; 4]; 4],
    columns: &[u32],
    images: &[[[u32; 4]; 4]],
) -> FiberScreen {
    screen_fibers_scalar(cursor, offsets, columns, images, 16)
}
struct FiberEnumeration {
    model: Option<u64>,
    complete: bool,
    logical: Logical,
    trace: u64,
}
fn enumerate_fibers(plan: &FiberPlan, cap: u64, mode: u8) -> FiberEnumeration {
    assert!(mode <= 2);
    assert!(plan.n <= 24 && plan.low.len() + plan.outside.len() == plan.n);
    let h = plan.outside.len();
    let size = 1usize << h.min(4);
    let high = h.saturating_sub(4);
    let mut cursor = GrayCursor::new(&plan.form);
    let offsets = plan.form.low_offsets();
    let mut columns = plan.columns.clone();
    let mut result = FiberEnumeration {
        model: None,
        complete: false,
        logical: Logical {
            fiber_selected: u64::from(plan.selected),
            fiber_dimension: plan.low.len() as u64,
            ..Logical::default()
        },
        trace: 0xcbf29ce484222325 ^ u64::from(plan.selected),
    };
    for step in 0..1u64 << high {
        if cap - result.logical.fiber_prefixes < size as u64 {
            return result;
        }
        cursor.advance(&plan.form, step);
        if step != 0 {
            let j = step.trailing_zeros() as usize;
            for i in 0..columns.len() {
                columns[i] ^= plan.cross_high[j][i];
            }
        }
        let screen = if mode == 2 && size == 16 {
            screen_fibers_native(&cursor, &offsets, &columns, &plan.column_offsets)
        } else {
            screen_fibers_scalar(&cursor, &offsets, &columns, &plan.column_offsets, size)
        };
        result.logical.fiber_prefixes += size as u64;
        result.logical.fiber_filter_rounds += screen.rounds;
        result.logical.fiber_batches += 1;
        let rejected = size as u64 - u64::from(screen.survivors.count_ones());
        result.logical.fiber_zero_rejected += rejected;
        result.logical.fiber_extensions_rejected += rejected << plan.low.len();
        result.trace = (result.trace ^ step).wrapping_mul(0x100000001b3);
        result.trace = (result.trace ^ u64::from(screen.survivors)).wrapping_mul(0x100000001b3);
        let mut survivors = screen.survivors;
        while survivors != 0 {
            let y = survivors.trailing_zeros() as usize;
            survivors &= survivors - 1;
            let mut a = [0u32; MAX];
            for i in 0..columns.len() {
                a[i] = columns[i] ^ plan.column_offsets[i][y / 4][y % 4];
            }
            let b = screen.values[y / 4][y % 4];
            let solved = if mode == 0 {
                fiber_rows(&a[..columns.len()], b, plan.equations)
            } else {
                fiber_columns(&a[..columns.len()], b)
            };
            result.logical.fiber_queries += 1;
            result.logical.fiber_rank_sum += u64::from(solved.rank);
            result.trace = (result.trace ^ y as u64).wrapping_mul(0x100000001b3);
            result.trace = (result.trace ^ u64::from(solved.rank)).wrapping_mul(0x100000001b3);
            result.trace = (result.trace ^ solved.model.map_or(0, |m| u64::from(m) + 1))
                .wrapping_mul(0x100000001b3);
            if let Some(inside) = solved.model {
                let outside_point = ((step ^ (step >> 1)) << 4) | y as u64;
                result.model = Some(plan.recover(outside_point, inside));
                result.complete = true;
                return result;
            }
            result.logical.fiber_linear_rejected += 1;
            result.logical.fiber_extensions_rejected += 1u64 << plan.low.len();
        }
    }
    result.complete = true;
    result
}
fn solve_fiber(system: &System, n: u8, limit: u64, mode: u8) -> Solved {
    if n > 24 {
        return solve_packed(system, n, limit);
    }
    let tick = Instant::now();
    let Some(original) = SyndromeForm::from_system(system, n) else {
        return solve_packed(system, n, limit);
    };
    let graph = fiber_graph(&original);
    let (selected, states) = maximum_independent(&graph, n as usize);
    let selection_ns = tick.elapsed().as_nanos();
    let tick = Instant::now();
    let plan = FiberPlan::new(&original, system.len(), selected);
    let fiber_setup_ns = tick.elapsed().as_nanos();
    let tick = Instant::now();
    let mut got = enumerate_fibers(&plan, if limit == 0 { 0 } else { 1 << 24 }, mode);
    let fiber_ns = tick.elapsed().as_nanos();
    got.logical.fiber_selection_states = states;
    Solved {
        outcome: if !got.complete {
            Outcome::Unknown("FIBER_CAP")
        } else {
            got.model.map_or(Outcome::Unsat, Outcome::Sat)
        },
        logical: got.logical,
        profile: Profile {
            selection_ns,
            fiber_setup_ns,
            fiber_ns,
            kernel_calls_by_active: vec![0; n as usize + 1],
            ..Profile::default()
        },
        trace: got.trace,
    }
}
