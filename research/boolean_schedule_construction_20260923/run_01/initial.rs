// One complete affine tail of the ORIGINAL coefficient row space, followed by
// one injective parametrization. No polynomial-product closure is claimed.
fn projection_slot(monomial: u64) -> usize {
    (monomial.wrapping_mul(0x9e3779b97f4a7c15).rotate_left(23) >> 58) as usize
}
fn certify_quadratic_row_rank(system: &System) -> (bool, u64, u64) {
    // A fixed LINEAR map from the quadratic coefficient space into F_2^64.
    // Full row rank after projection certifies full row rank before projection.
    // A hash collision can only make the certificate fail; failure falls back.
    let mut pivots = [0u64; 64];
    let mut rows = 0;
    let mut xors = 0;
    for poly in system {
        rows += 1;
        let mut row = 0;
        for &m in poly {
            if m.count_ones() == 2 {
                row ^= 1u64 << projection_slot(m);
            }
        }
        loop {
            if row == 0 {
                return (false, rows, xors);
            }
            let p = row.trailing_zeros() as usize;
            if pivots[p] == 0 {
                pivots[p] = row;
                break;
            }
            row ^= pivots[p];
            xors += 1;
        }
    }
    (true, rows, xors)
}
fn initial_affine_rows<B: BasisBackend>(backend: B, system: &System, n: u8) -> Vec<u64> {
    let state = backend.reduce(backend.build(system));
    let mut affine = Vec::new();
    for r in 0..backend.row_count(&state) {
        let mut row = 0;
        let mut nonlinear = false;
        backend.terms(&state, r, |m| {
            if m.count_ones() > 1 {
                nonlinear = true;
            } else {
                row ^= if m == 0 { 1u64 << n } else { m };
            }
        });
        if !nonlinear {
            affine.push(row);
        }
    }
    algebra::affine_rref(n, affine.into_iter())
}
struct InitialMap {
    images: RecoveryMap,
    free: u8,
    rank: usize,
}
impl InitialMap {
    fn from_rows(rows: &[u64], n: u8) -> Option<Self> {
        let constant = 1u64 << n;
        if rows.contains(&constant) {
            return None;
        }
        let mut images = identity_affine(constant - 1, n);
        let mut pivots = 0;
        for &row in rows {
            let variables = row & (constant - 1);
            if variables != 0 {
                let p = variables & variables.wrapping_neg();
                images[p.trailing_zeros() as usize] = row ^ p;
                pivots |= p;
            }
        }
        let free = n - pivots.count_ones() as u8;
        let mut labels = [0usize; MAX];
        let mut next = 0;
        for j in 0..n as usize {
            if pivots & (1u64 << j) == 0 {
                labels[j] = next;
                next += 1;
            }
        }
        for image in &mut images[..n as usize] {
            assert_eq!(*image & pivots, 0);
            let mut result = if *image & constant != 0 {
                1u64 << free
            } else {
                0
            };
            let mut vars = *image & (constant - 1);
            while vars != 0 {
                result ^= 1u64 << labels[vars.trailing_zeros() as usize];
                vars &= vars - 1;
            }
            *image = result;
        }
        Some(Self {
            images,
            free,
            rank: pivots.count_ones() as usize,
        })
    }
    fn recover(&self, point: u64, n: u8) -> u64 {
        let mut result = 0;
        for j in 0..n as usize {
            let bit = (self.images[j] & point).count_ones() % 2
                ^ ((self.images[j] >> self.free) & 1) as u32;
            result |= u64::from(bit) << j;
        }
        result
    }
}
fn add_affine_syndrome(form: &mut SyndromeForm, image: u64, coefficient: u32) {
    if image & (1u64 << form.n) != 0 {
        form.constant ^= coefficient;
    }
    let mut variables = image & ((1u64 << form.n) - 1);
    while variables != 0 {
        form.linear[variables.trailing_zeros() as usize] ^= coefficient;
        variables &= variables - 1;
    }
}
fn transform_syndrome(original: &SyndromeForm, map: &InitialMap) -> SyndromeForm {
    let mut out = SyndromeForm::new(map.free as usize);
    out.constant = original.constant;
    let constant = 1u64 << map.free;
    for i in 0..original.n {
        add_affine_syndrome(&mut out, map.images[i], original.linear[i]);
        for j in i + 1..original.n {
            let coefficient = original.quadratic[i][j];
            if coefficient == 0 {
                continue;
            }
            let (a, b) = (map.images[i], map.images[j]);
            // Constant-linear cross terms plus Boolean y_i*y_i=y_i.
            if a & constant != 0 {
                add_affine_syndrome(&mut out, b, coefficient);
            }
            if b & constant != 0 {
                add_affine_syndrome(&mut out, a & (constant - 1), coefficient);
            }
            let mut av = a & (constant - 1);
            while av != 0 {
                let x = av.trailing_zeros() as usize;
                av &= av - 1;
                let mut bv = b & (constant - 1);
                while bv != 0 {
                    let y = bv.trailing_zeros() as usize;
                    bv &= bv - 1;
                    if x == y {
                        out.linear[x] ^= coefficient;
                    } else {
                        out.quadratic[x][y] ^= coefficient;
                        out.quadratic[y][x] ^= coefficient;
                    }
                }
            }
        }
    }
    out
}
fn solve_initial(system: &System, n: u8, limit: u64, optimized: bool) -> Solved {
    if n > 24 {
        return solve_packed(system, n, limit);
    }
    let tick = Instant::now();
    let Some(original) = SyndromeForm::from_system(system, n) else {
        return solve_packed(system, n, limit);
    };
    if limit == 0 {
        return solve_gray_delta(system, n, limit, optimized);
    }
    let (certified, projection_rows, projection_xors) = certify_quadratic_row_rank(system);
    let rows = if certified {
        vec![]
    } else {
        let mut canonical = system.clone();
        for p in &mut canonical {
            algebra::canonical(p);
        }
        if optimized {
            initial_affine_rows(WideBasis::tail(n), &canonical, n)
        } else {
            initial_affine_rows(TailList, &canonical, n)
        }
    };
    let calls = if certified { 1 } else { 2 };
    let mut profile = Profile {
        kernel_ns: tick.elapsed().as_nanos(),
        kernel_calls_by_active: vec![0; n as usize + 1],
        ..Profile::default()
    };
    profile.kernel_calls_by_active[n as usize] = calls;
    let mut logical = Logical {
        kernel_calls: calls,
        source_rows: projection_rows + if certified { 0 } else { system.len() as u64 },
        source_columns: 64
            + if certified {
                0
            } else {
                (usize::from(n) * usize::from(n.saturating_sub(1)) / 2 + usize::from(n) + 1) as u64
            },
        projection_rows,
        projection_xors,
        projection_certified: u64::from(certified),
        specialized_terms: system.iter().map(|p| p.len() as u64).sum(),
        derived_rows: rows.len() as u64,
        ..Logical::default()
    };
    let tick = Instant::now();
    let Some(map) = InitialMap::from_rows(&rows, n) else {
        profile.substitution_ns = tick.elapsed().as_nanos();
        return Solved {
            outcome: Outcome::Unsat,
            logical,
            profile,
            trace: 0xaff1bad,
        };
    };
    logical.affine_eliminated = map.rank as u64;
    let form = if map.rank == 0 {
        original
    } else if optimized {
        transform_syndrome(&original, &map)
    } else {
        let transformed: System = system
            .iter()
            .map(|p| substitute_polynomial(p, &map.images, map.free))
            .collect();
        SyndromeForm::from_system(&transformed, map.free).unwrap()
    };
    profile.substitution_ns = tick.elapsed().as_nanos();
    let tick = Instant::now();
    let got = if optimized {
        enumerate_delta::<NativeDeltaBlock>(&form, 1 << 24)
    } else {
        enumerate_delta::<ScalarDeltaBlock>(&form, 1 << 24)
    };
    profile.enumeration_ns = tick.elapsed().as_nanos();
    logical.enumeration_points = got.points;
    logical.enumeration_batches = got.batches;
    logical.enumeration_leaves = 1;
    let mut trace = got.checksum;
    for image in &map.images[..n as usize] {
        trace = (trace ^ image).wrapping_mul(0x100000001b3);
    }
    let outcome = if !got.complete {
        Outcome::Unknown("ENUM_CAP")
    } else {
        got.model
            .map_or(Outcome::Unsat, |point| Outcome::Sat(map.recover(point, n)))
    };
    Solved {
        outcome,
        logical,
        profile,
        trace,
    }
}
