// Transport all sixteen equation-syndrome values through high Gray steps.
// The changing difference is affine in the four low coordinates. When high
// coordinate j changes, lower scheduled differences change by a uniform word.
trait DeltaBlock: Copy {
    fn from_values(values: &[[u32; 4]; 4]) -> Self;
    fn values(&self) -> [[u32; 4]; 4];
    fn xor_block(&mut self, other: &Self);
    fn xor_uniform(&mut self, word: u32);
    fn inspect(&self) -> (u32, Option<u8>);
}
#[derive(Clone, Copy)]
struct ScalarDeltaBlock([[u32; 4]; 4]);
fn inspect_delta_values(values: &[[u32; 4]; 4]) -> (u32, Option<u8>) {
    let mut sum = 0u32;
    let mut found = None;
    for (y, &value) in values.iter().flatten().enumerate() {
        sum = sum.wrapping_add(value);
        if value == 0 && found.is_none() {
            found = Some(y as u8);
        }
    }
    (sum, found)
}
impl DeltaBlock for ScalarDeltaBlock {
    fn from_values(values: &[[u32; 4]; 4]) -> Self {
        Self(*values)
    }
    fn values(&self) -> [[u32; 4]; 4] {
        self.0
    }
    fn xor_block(&mut self, other: &Self) {
        for i in 0..4 {
            for j in 0..4 {
                self.0[i][j] ^= other.0[i][j];
            }
        }
    }
    fn xor_uniform(&mut self, word: u32) {
        for value in self.0.iter_mut().flatten() {
            *value ^= word;
        }
    }
    fn inspect(&self) -> (u32, Option<u8>) {
        inspect_delta_values(&self.0)
    }
}

#[cfg(target_arch = "aarch64")]
#[derive(Clone, Copy)]
struct NativeDeltaBlock([std::arch::aarch64::uint32x4_t; 4]);
#[cfg(target_arch = "aarch64")]
impl DeltaBlock for NativeDeltaBlock {
    fn from_values(values: &[[u32; 4]; 4]) -> Self {
        // Every load addresses a complete four-u32 array. NEON is baseline.
        unsafe {
            Self(std::array::from_fn(|i| {
                std::arch::aarch64::vld1q_u32(values[i].as_ptr())
            }))
        }
    }
    fn values(&self) -> [[u32; 4]; 4] {
        let mut out = [[0; 4]; 4];
        for i in 0..4 {
            unsafe {
                std::arch::aarch64::vst1q_u32(out[i].as_mut_ptr(), self.0[i]);
            }
        }
        out
    }
    fn xor_block(&mut self, other: &Self) {
        for i in 0..4 {
            self.0[i] = unsafe { std::arch::aarch64::veorq_u32(self.0[i], other.0[i]) };
        }
    }
    fn xor_uniform(&mut self, word: u32) {
        unsafe {
            use std::arch::aarch64::*;
            let value = vdupq_n_u32(word);
            for lane in &mut self.0 {
                *lane = veorq_u32(*lane, value);
            }
        }
    }
    fn inspect(&self) -> (u32, Option<u8>) {
        unsafe {
            use std::arch::aarch64::*;
            let [a, b, c, d] = self.0;
            let sum = vaddvq_u32(vaddq_u32(vaddq_u32(a, b), vaddq_u32(c, d)));
            if vminvq_u32(vminq_u32(vminq_u32(a, b), vminq_u32(c, d))) == 0 {
                let checked = inspect_delta_values(&self.values());
                assert_eq!(checked.0, sum);
                checked
            } else {
                (sum, None)
            }
        }
    }
}
#[cfg(target_arch = "x86_64")]
#[derive(Clone, Copy)]
struct NativeDeltaBlock([std::arch::x86_64::__m128i; 4]);
#[cfg(target_arch = "x86_64")]
impl DeltaBlock for NativeDeltaBlock {
    fn from_values(values: &[[u32; 4]; 4]) -> Self {
        unsafe {
            Self(std::array::from_fn(|i| {
                std::arch::x86_64::_mm_loadu_si128(values[i].as_ptr().cast())
            }))
        }
    }
    fn values(&self) -> [[u32; 4]; 4] {
        let mut out = [[0; 4]; 4];
        for i in 0..4 {
            unsafe {
                std::arch::x86_64::_mm_storeu_si128(out[i].as_mut_ptr().cast(), self.0[i]);
            }
        }
        out
    }
    fn xor_block(&mut self, other: &Self) {
        for i in 0..4 {
            self.0[i] = unsafe { std::arch::x86_64::_mm_xor_si128(self.0[i], other.0[i]) };
        }
    }
    fn xor_uniform(&mut self, word: u32) {
        unsafe {
            use std::arch::x86_64::*;
            let value = _mm_set1_epi32(word as i32);
            for lane in &mut self.0 {
                *lane = _mm_xor_si128(*lane, value);
            }
        }
    }
    fn inspect(&self) -> (u32, Option<u8>) {
        unsafe {
            use std::arch::x86_64::*;
            let [a, b, c, d] = self.0;
            let mut sum = _mm_add_epi32(_mm_add_epi32(a, b), _mm_add_epi32(c, d));
            sum = _mm_add_epi32(sum, _mm_shuffle_epi32(sum, 0x4e));
            sum = _mm_add_epi32(sum, _mm_shuffle_epi32(sum, 0xb1));
            let total = _mm_cvtsi128_si32(sum) as u32;
            let z = _mm_setzero_si128();
            let equal = _mm_or_si128(
                _mm_or_si128(_mm_cmpeq_epi32(a, z), _mm_cmpeq_epi32(b, z)),
                _mm_or_si128(_mm_cmpeq_epi32(c, z), _mm_cmpeq_epi32(d, z)),
            );
            if _mm_movemask_epi8(equal) != 0 {
                let checked = inspect_delta_values(&self.values());
                assert_eq!(checked.0, total);
                checked
            } else {
                (total, None)
            }
        }
    }
}
#[cfg(not(any(target_arch = "aarch64", target_arch = "x86_64")))]
type NativeDeltaBlock = ScalarDeltaBlock;

struct DeltaCursor<B: DeltaBlock> {
    block: B,
    differences: [B; MAX],
}
impl<B: DeltaBlock> DeltaCursor<B> {
    fn new(form: &SyndromeForm) -> Self {
        let zero = B::from_values(&[[0; 4]; 4]);
        let mut values = [[0; 4]; 4];
        for y in 0..16 {
            values[y / 4][y % 4] = form.value(y as u64);
        }
        let mut differences = [zero; MAX];
        for j in 0..form.n - 4 {
            let c = form.linear[j + 4]
                ^ if j == 0 {
                    0
                } else {
                    form.quadratic[j + 4][j + 3]
                };
            let mut delta = [[c; 4]; 4];
            for y in 0..16 {
                for i in 0..4 {
                    if y & (1 << i) != 0 {
                        delta[y / 4][y % 4] ^= form.quadratic[i][j + 4];
                    }
                }
            }
            differences[j] = B::from_values(&delta);
        }
        Self {
            block: B::from_values(&values),
            differences,
        }
    }
    fn advance(&mut self, form: &SyndromeForm, step: u64) {
        if step == 0 {
            return;
        }
        let j = step.trailing_zeros() as usize;
        self.block.xor_block(&self.differences[j]);
        for i in 0..j {
            self.differences[i].xor_uniform(form.quadratic[i + 4][j + 4]);
        }
    }
}
fn enumerate_delta<B: DeltaBlock>(form: &SyndromeForm, cap: u64) -> Enumeration {
    if form.n < 4 {
        return enumerate_syndromes::<false>(form, cap);
    }
    let mut result = Enumeration {
        model: None,
        complete: false,
        points: 0,
        batches: 0,
        checksum: 0xcbf29ce484222325,
    };
    let mut cursor = DeltaCursor::<B>::new(form);
    for step in 0..1u64 << (form.n - 4) {
        if cap - result.points < 16 {
            return result;
        }
        cursor.advance(form, step);
        let (sum, found) = cursor.block.inspect();
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
fn solve_gray_delta(system: &System, n: u8, limit: u64, simd: bool) -> Solved {
    if n > 24 {
        return solve_packed(system, n, limit);
    }
    let tick = Instant::now();
    let Some(form) = SyndromeForm::from_system(system, n) else {
        return solve_packed(system, n, limit);
    };
    let cap = if limit == 0 { 0 } else { 1 << 24 };
    let got = if simd {
        enumerate_delta::<NativeDeltaBlock>(&form, cap)
    } else {
        enumerate_delta::<ScalarDeltaBlock>(&form, cap)
    };
    Solved {
        outcome: if !got.complete {
            Outcome::Unknown("ENUM_CAP")
        } else {
            got.model.map_or(Outcome::Unsat, Outcome::Sat)
        },
        logical: Logical {
            enumeration_points: got.points,
            enumeration_batches: got.batches,
            enumeration_leaves: 1,
            ..Logical::default()
        },
        profile: Profile {
            enumeration_ns: tick.elapsed().as_nanos(),
            kernel_calls_by_active: vec![0; n as usize + 1],
            ..Profile::default()
        },
        trace: got.checksum,
    }
}
