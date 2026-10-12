// Each bit is one equation value. Gray differences update the high-variable
// assignment; four low variables are evaluated as a complete 16-point block.
struct SyndromeForm {
    n: usize,
    constant: u32,
    linear: [u32; MAX],
    quadratic: [[u32; MAX]; MAX],
}
impl SyndromeForm {
    fn new(n: usize) -> Self {
        Self {
            n,
            constant: 0,
            linear: [0; MAX],
            quadratic: [[0; MAX]; MAX],
        }
    }
    fn add_term(&mut self, equation: usize, monomial: u64) {
        let bit = 1u32 << equation;
        match monomial.count_ones() {
            0 => self.constant ^= bit,
            1 => self.linear[monomial.trailing_zeros() as usize] ^= bit,
            2 => {
                let a = monomial.trailing_zeros() as usize;
                let b = (monomial & (monomial - 1)).trailing_zeros() as usize;
                self.quadratic[a][b] ^= bit;
                self.quadratic[b][a] ^= bit;
            }
            _ => unreachable!(),
        }
    }
    fn from_system(system: &System, n: u8) -> Option<Self> {
        if n as usize > MAX
            || system.len() > 32
            || system
                .iter()
                .any(|p| p.iter().any(|&m| m >= 1u64 << n || m.count_ones() > 2))
        {
            return None;
        }
        let mut form = Self::new(n as usize);
        for (e, p) in system.iter().enumerate() {
            for &m in p {
                form.add_term(e, m);
            }
        }
        Some(form)
    }
    fn value(&self, point: u64) -> u32 {
        let mut value = self.constant;
        for i in 0..self.n {
            if point & (1u64 << i) != 0 {
                value ^= self.linear[i];
                for j in i + 1..self.n {
                    if point & (1u64 << j) != 0 {
                        value ^= self.quadratic[i][j];
                    }
                }
            }
        }
        value
    }
    fn low_offsets(&self) -> [[u32; 4]; 4] {
        let mut offsets = [[0; 4]; 4];
        for y in 0..16 {
            for i in 0..4 {
                for j in i + 1..4 {
                    if y & (1 << i) != 0 && y & (1 << j) != 0 {
                        offsets[y / 4][y % 4] ^= self.quadratic[i][j];
                    }
                }
            }
        }
        offsets
    }
}
fn scalar_syndrome_block(
    constant: u32,
    linear: &[u32; 4],
    offsets: &[[u32; 4]; 4],
) -> (u32, Option<u8>) {
    let mut sum = 0u32;
    let mut found = None;
    for y in 0..16 {
        let mut value = constant ^ offsets[y / 4][y % 4];
        for i in 0..4 {
            if y & (1 << i) != 0 {
                value ^= linear[i];
            }
        }
        sum = sum.wrapping_add(value);
        if value == 0 && found.is_none() {
            found = Some(y as u8);
        }
    }
    (sum, found)
}

#[cfg(target_arch = "aarch64")]
fn simd_syndrome_block(
    constant: u32,
    linear: &[u32; 4],
    offsets: &[[u32; 4]; 4],
) -> (u32, Option<u8>) {
    // AArch64 guarantees NEON. Every load addresses one complete four-u32 array.
    unsafe {
        use std::arch::aarch64::*;
        let base = [
            constant,
            constant ^ linear[0],
            constant ^ linear[1],
            constant ^ linear[0] ^ linear[1],
        ];
        let v = vld1q_u32(base.as_ptr());
        let x2 = vdupq_n_u32(linear[2]);
        let x3 = vdupq_n_u32(linear[3]);
        let a = veorq_u32(v, vld1q_u32(offsets[0].as_ptr()));
        let b = veorq_u32(veorq_u32(v, x2), vld1q_u32(offsets[1].as_ptr()));
        let c = veorq_u32(veorq_u32(v, x3), vld1q_u32(offsets[2].as_ptr()));
        let d = veorq_u32(
            veorq_u32(veorq_u32(v, x2), x3),
            vld1q_u32(offsets[3].as_ptr()),
        );
        let sum = vaddvq_u32(vaddq_u32(vaddq_u32(a, b), vaddq_u32(c, d)));
        let any = vminvq_u32(vminq_u32(vminq_u32(a, b), vminq_u32(c, d))) == 0;
        if any {
            let checked = scalar_syndrome_block(constant, linear, offsets);
            assert_eq!(checked.0, sum);
            checked
        } else {
            (sum, None)
        }
    }
}
#[cfg(target_arch = "x86_64")]
fn simd_syndrome_block(
    constant: u32,
    linear: &[u32; 4],
    offsets: &[[u32; 4]; 4],
) -> (u32, Option<u8>) {
    // SSE2 is baseline on x86_64; unaligned loads read complete four-u32 arrays.
    unsafe {
        use std::arch::x86_64::*;
        let base = [
            constant,
            constant ^ linear[0],
            constant ^ linear[1],
            constant ^ linear[0] ^ linear[1],
        ];
        let v = _mm_loadu_si128(base.as_ptr().cast());
        let x2 = _mm_set1_epi32(linear[2] as i32);
        let x3 = _mm_set1_epi32(linear[3] as i32);
        let a = _mm_xor_si128(v, _mm_loadu_si128(offsets[0].as_ptr().cast()));
        let b = _mm_xor_si128(
            _mm_xor_si128(v, x2),
            _mm_loadu_si128(offsets[1].as_ptr().cast()),
        );
        let c = _mm_xor_si128(
            _mm_xor_si128(v, x3),
            _mm_loadu_si128(offsets[2].as_ptr().cast()),
        );
        let d = _mm_xor_si128(
            _mm_xor_si128(_mm_xor_si128(v, x2), x3),
            _mm_loadu_si128(offsets[3].as_ptr().cast()),
        );
        let mut sum = _mm_add_epi32(_mm_add_epi32(a, b), _mm_add_epi32(c, d));
        sum = _mm_add_epi32(sum, _mm_shuffle_epi32(sum, 0x4e));
        sum = _mm_add_epi32(sum, _mm_shuffle_epi32(sum, 0xb1));
        let total = _mm_cvtsi128_si32(sum) as u32;
        let zero = _mm_setzero_si128();
        let matches = _mm_or_si128(
            _mm_or_si128(_mm_cmpeq_epi32(a, zero), _mm_cmpeq_epi32(b, zero)),
            _mm_or_si128(_mm_cmpeq_epi32(c, zero), _mm_cmpeq_epi32(d, zero)),
        );
        if _mm_movemask_epi8(matches) != 0 {
            let checked = scalar_syndrome_block(constant, linear, offsets);
            assert_eq!(checked.0, total);
            checked
        } else {
            (total, None)
        }
    }
}
#[cfg(not(any(target_arch = "aarch64", target_arch = "x86_64")))]
fn simd_syndrome_block(
    constant: u32,
    linear: &[u32; 4],
    offsets: &[[u32; 4]; 4],
) -> (u32, Option<u8>) {
    scalar_syndrome_block(constant, linear, offsets)
}

struct GrayCursor {
    constant: u32,
    linear: [u32; 4],
    differences: [u32; MAX],
}
impl GrayCursor {
    fn new(form: &SyndromeForm) -> Self {
        let mut differences = [0u32; MAX];
        for j in 0..form.n.saturating_sub(4) {
            differences[j] = form.linear[j + 4]
                ^ if j == 0 {
                    0
                } else {
                    form.quadratic[j + 4][j + 3]
                };
        }
        Self {
            constant: form.constant,
            linear: [
                form.linear[0],
                form.linear[1],
                form.linear[2],
                form.linear[3],
            ],
            differences,
        }
    }
    fn advance(&mut self, form: &SyndromeForm, step: u64) {
        if step == 0 {
            return;
        }
        let j = step.trailing_zeros() as usize;
        self.constant ^= self.differences[j];
        for i in 0..j {
            self.differences[i] ^= form.quadratic[i + 4][j + 4];
        }
        for i in 0..4 {
            self.linear[i] ^= form.quadratic[i][j + 4];
        }
    }
}
struct Enumeration {
    model: Option<u64>,
    complete: bool,
    points: u64,
    batches: u64,
    checksum: u64,
}
fn enumerate_syndromes<const SIMD: bool>(form: &SyndromeForm, point_cap: u64) -> Enumeration {
    let mut result = Enumeration {
        model: None,
        complete: false,
        points: 0,
        batches: 0,
        checksum: 0xcbf29ce484222325,
    };
    if form.n < 4 {
        for point in 0..1u64 << form.n {
            if result.points == point_cap {
                return result;
            }
            let value = form.value(point);
            result.points += 1;
            result.batches += 1;
            result.checksum = (result.checksum ^ u64::from(value)).wrapping_mul(0x100000001b3);
            if value == 0 {
                result.model = Some(point);
                result.complete = true;
                return result;
            }
        }
        result.complete = true;
        return result;
    }
    let high = form.n - 4;
    let mut cursor = GrayCursor::new(form);
    let offsets = form.low_offsets();
    for step in 0..1u64 << high {
        if point_cap - result.points < 16 {
            return result;
        }
        cursor.advance(form, step);
        let (sum, found) = if SIMD {
            simd_syndrome_block(cursor.constant, &cursor.linear, &offsets)
        } else {
            scalar_syndrome_block(cursor.constant, &cursor.linear, &offsets)
        };
        result.points += 16;
        result.batches += 1;
        result.checksum = (result.checksum ^ step).wrapping_mul(0x100000001b3);
        result.checksum = (result.checksum ^ u64::from(sum)).wrapping_mul(0x100000001b3);
        if let Some(low) = found {
            result.model = Some(((step ^ (step >> 1)) << 4) | u64::from(low));
            result.complete = true;
            return result;
        }
    }
    result.complete = true;
    result
}
fn solve_gray(system: &System, n: u8, limit: u64, simd: bool) -> Solved {
    if n > 24 {
        return solve_packed(system, n, limit);
    }
    let start = Instant::now();
    let Some(form) = SyndromeForm::from_system(system, n) else {
        return solve_packed(system, n, limit);
    };
    let cap = if limit == 0 { 0 } else { 1u64 << 24 };
    let got = if simd {
        enumerate_syndromes::<true>(&form, cap)
    } else {
        enumerate_syndromes::<false>(&form, cap)
    };
    let enumeration_ns = start.elapsed().as_nanos();
    let outcome = if !got.complete {
        Outcome::Unknown("ENUM_CAP")
    } else {
        got.model.map_or(Outcome::Unsat, Outcome::Sat)
    };
    Solved {
        outcome,
        logical: Logical {
            enumeration_points: got.points,
            enumeration_batches: got.batches,
            enumeration_leaves: 1,
            ..Logical::default()
        },
        profile: Profile {
            enumeration_ns,
            kernel_calls_by_active: vec![0; n as usize + 1],
            ..Profile::default()
        },
        trace: got.checksum,
    }
}
