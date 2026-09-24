// Exact projection onto the first eight equation bits, with full checks on
// survivors. Universal masks contain no fixture-dependent preprocessing.
const fn byte_low_masks() -> [[[u8; 16]; 6]; 4] {
    let mut out = [[[0; 16]; 6]; 4];
    let mut g = 0;
    while g < 4 {
        let mut i = 0;
        while i < 6 {
            let mut lane = 0;
            while lane < 16 {
                out[g][i][lane] = if ((g * 16 + lane) & (1 << i)) != 0 {
                    255
                } else {
                    0
                };
                lane += 1;
            }
            i += 1;
        }
        g += 1;
    }
    out
}
const BYTE_LOW_MASKS: [[[u8; 16]; 6]; 4] = byte_low_masks();
trait ByteBlock64: Copy {
    fn from_values(values: &[u8; 64]) -> Self;
    fn affine(constant: u8, linear: &[u8; 6]) -> Self;
    fn values(&self) -> [u8; 64];
    fn xor_uniform(&mut self, value: u8);
    fn xor_block(&mut self, other: &Self);
    fn zero_masks(&self) -> [u16; 4];
}
#[derive(Clone, Copy)]
struct ScalarByte64([u8; 64]);
impl ByteBlock64 for ScalarByte64 {
    fn from_values(values: &[u8; 64]) -> Self {
        Self(*values)
    }
    fn affine(constant: u8, linear: &[u8; 6]) -> Self {
        let mut values = [constant; 64];
        for y in 1usize..64 {
            let bit = y & y.wrapping_neg();
            values[y] = values[y ^ bit] ^ linear[bit.trailing_zeros() as usize];
        }
        Self(values)
    }
    fn values(&self) -> [u8; 64] {
        self.0
    }
    fn xor_uniform(&mut self, value: u8) {
        for v in &mut self.0 {
            *v ^= value;
        }
    }
    fn xor_block(&mut self, other: &Self) {
        for i in 0..64 {
            self.0[i] ^= other.0[i];
        }
    }
    fn zero_masks(&self) -> [u16; 4] {
        let mut masks = [0; 4];
        for y in 0..64 {
            if self.0[y] == 0 {
                masks[y / 16] |= 1 << (y % 16);
            }
        }
        masks
    }
}
#[cfg(target_arch = "aarch64")]
#[derive(Clone, Copy)]
struct NativeByte64([std::arch::aarch64::uint8x16_t; 4]);
#[cfg(target_arch = "aarch64")]
impl ByteBlock64 for NativeByte64 {
    fn from_values(values: &[u8; 64]) -> Self {
        unsafe {
            Self(std::array::from_fn(|g| {
                std::arch::aarch64::vld1q_u8(values[g * 16..].as_ptr())
            }))
        }
    }
    fn affine(constant: u8, linear: &[u8; 6]) -> Self {
        unsafe {
            use std::arch::aarch64::*;
            Self(std::array::from_fn(|g| {
                let mut v = vdupq_n_u8(constant);
                for i in 0..6 {
                    v = veorq_u8(
                        v,
                        vandq_u8(
                            vdupq_n_u8(linear[i]),
                            vld1q_u8(BYTE_LOW_MASKS[g][i].as_ptr()),
                        ),
                    );
                }
                v
            }))
        }
    }
    fn values(&self) -> [u8; 64] {
        let mut out = [0; 64];
        for g in 0..4 {
            unsafe {
                std::arch::aarch64::vst1q_u8(out[g * 16..].as_mut_ptr(), self.0[g]);
            }
        }
        out
    }
    fn xor_uniform(&mut self, value: u8) {
        unsafe {
            let v = std::arch::aarch64::vdupq_n_u8(value);
            for x in &mut self.0 {
                *x = std::arch::aarch64::veorq_u8(*x, v);
            }
        }
    }
    fn xor_block(&mut self, other: &Self) {
        for g in 0..4 {
            self.0[g] = unsafe { std::arch::aarch64::veorq_u8(self.0[g], other.0[g]) };
        }
    }
    fn zero_masks(&self) -> [u16; 4] {
        unsafe {
            use std::arch::aarch64::*;
            let [a, b, c, d] = self.0;
            let mut masks = [0; 4];
            if vminvq_u8(vminq_u8(vminq_u8(a, b), vminq_u8(c, d))) != 0 {
                return masks;
            }
            let weights = [1u8, 2, 4, 8, 16, 32, 64, 128, 1, 2, 4, 8, 16, 32, 64, 128];
            let weights = vld1q_u8(weights.as_ptr());
            for g in 0..4 {
                if vminvq_u8(self.0[g]) == 0 {
                    let selected = vandq_u8(vceqq_u8(self.0[g], vdupq_n_u8(0)), weights);
                    masks[g] = u16::from(vaddv_u8(vget_low_u8(selected)))
                        | (u16::from(vaddv_u8(vget_high_u8(selected))) << 8);
                }
            }
            masks
        }
    }
}
#[cfg(target_arch = "x86_64")]
#[derive(Clone, Copy)]
struct NativeByte64([std::arch::x86_64::__m128i; 4]);
#[cfg(target_arch = "x86_64")]
impl ByteBlock64 for NativeByte64 {
    fn from_values(values: &[u8; 64]) -> Self {
        unsafe {
            Self(std::array::from_fn(|g| {
                std::arch::x86_64::_mm_loadu_si128(values[g * 16..].as_ptr().cast())
            }))
        }
    }
    fn affine(constant: u8, linear: &[u8; 6]) -> Self {
        unsafe {
            use std::arch::x86_64::*;
            Self(std::array::from_fn(|g| {
                let mut v = _mm_set1_epi8(constant as i8);
                for i in 0..6 {
                    v = _mm_xor_si128(
                        v,
                        _mm_and_si128(
                            _mm_set1_epi8(linear[i] as i8),
                            _mm_loadu_si128(BYTE_LOW_MASKS[g][i].as_ptr().cast()),
                        ),
                    );
                }
                v
            }))
        }
    }
    fn values(&self) -> [u8; 64] {
        let mut out = [0; 64];
        for g in 0..4 {
            unsafe {
                std::arch::x86_64::_mm_storeu_si128(out[g * 16..].as_mut_ptr().cast(), self.0[g]);
            }
        }
        out
    }
    fn xor_uniform(&mut self, value: u8) {
        unsafe {
            let v = std::arch::x86_64::_mm_set1_epi8(value as i8);
            for x in &mut self.0 {
                *x = std::arch::x86_64::_mm_xor_si128(*x, v);
            }
        }
    }
    fn xor_block(&mut self, other: &Self) {
        for g in 0..4 {
            self.0[g] = unsafe { std::arch::x86_64::_mm_xor_si128(self.0[g], other.0[g]) };
        }
    }
    fn zero_masks(&self) -> [u16; 4] {
        unsafe {
            use std::arch::x86_64::*;
            let zero = _mm_setzero_si128();
            let tests = self.0.map(|v| _mm_cmpeq_epi8(v, zero));
            if _mm_movemask_epi8(_mm_or_si128(
                _mm_or_si128(tests[0], tests[1]),
                _mm_or_si128(tests[2], tests[3]),
            )) == 0
            {
                return [0; 4];
            }
            tests.map(|t| _mm_movemask_epi8(t) as u16)
        }
    }
}
#[cfg(not(any(target_arch = "aarch64", target_arch = "x86_64")))]
type NativeByte64 = ScalarByte64;

fn low_quadratic_offsets64(form: &SyndromeForm) -> [u32; 64] {
    let mut values = [0; 64];
    for y in 1usize..64 {
        let bit = y & y.wrapping_neg();
        let i = bit.trailing_zeros() as usize;
        let rest = y ^ bit;
        let mut value = values[rest];
        let mut bits = rest;
        while bits != 0 {
            let j = bits.trailing_zeros() as usize;
            bits &= bits - 1;
            value ^= form.quadratic[i][j];
        }
        values[y] = value;
    }
    values
}
struct SixCursor {
    constant: u32,
    linear: [u32; 6],
    differences: [u32; MAX],
    cross: Vec<[u32; 6]>,
}
impl SixCursor {
    fn new(form: &SyndromeForm) -> Self {
        let mut differences = [0; MAX];
        for j in 0..form.n - 6 {
            differences[j] = form.linear[j + 6]
                ^ if j == 0 {
                    0
                } else {
                    form.quadratic[j + 6][j + 5]
                };
        }
        Self {
            constant: form.constant,
            linear: std::array::from_fn(|i| form.linear[i]),
            differences,
            cross: (6..form.n)
                .map(|j| std::array::from_fn(|i| form.quadratic[i][j]))
                .collect(),
        }
    }
    #[inline(always)]
    fn advance<const LINEAR: bool>(
        &mut self,
        form: &SyndromeForm,
        step: u64,
    ) -> Option<(usize, u32)> {
        if step == 0 {
            return None;
        }
        let j = step.trailing_zeros() as usize;
        let delta = self.differences[j];
        self.constant ^= delta;
        for i in 0..j {
            self.differences[i] ^= form.quadratic[i + 6][j + 6];
        }
        if LINEAR {
            for i in 0..6 {
                self.linear[i] ^= self.cross[j][i];
            }
        }
        Some((j, delta))
    }
    fn value(&self, offsets: &[u32; 64], low: usize) -> u32 {
        let mut value = self.constant ^ offsets[low];
        let mut bits = low;
        while bits != 0 {
            let i = bits.trailing_zeros() as usize;
            bits &= bits - 1;
            value ^= self.linear[i];
        }
        value
    }
}
fn legacy_group64(step: u64, r: usize) -> usize {
    (r ^ (r >> 1)) ^ (((step as usize) & 1) << 1)
}
struct Filtered {
    model: Option<u64>,
    complete: bool,
    logical: Logical,
    trace: u64,
}
fn enumerate_bytes<B: ByteBlock64>(form: &SyndromeForm, cap: u64) -> Filtered {
    let mut result = Filtered {
        model: None,
        complete: false,
        logical: Logical {
            screen_calls: 1,
            ..Logical::default()
        },
        trace: 0xcbf29ce484222325,
    };
    if form.n < 6 {
        let got = enumerate_quiet(form, cap);
        result.model = got.model;
        result.complete = got.complete;
        result.logical.screen_fallback_points = got.points;
        result.logical.screen_fallback_batches = got.batches;
        result.trace ^= got.points ^ got.model.unwrap_or(u64::MAX);
        return result;
    }
    let offsets = low_quadratic_offsets64(form);
    let mut cursor = SixCursor::new(form);
    let mut block = B::affine(form.constant as u8, &cursor.linear.map(|v| v as u8));
    block.xor_block(&B::from_values(&offsets.map(|v| v as u8)));
    let cross: Vec<B> = cursor
        .cross
        .iter()
        .map(|a| B::affine(0, &a.map(|v| v as u8)))
        .collect();
    for step in 0..1u64 << (form.n - 6) {
        if cap - result.logical.screen_points < 64 {
            return result;
        }
        if let Some((j, d)) = cursor.advance::<true>(form, step) {
            block.xor_uniform(d as u8);
            block.xor_block(&cross[j]);
        }
        let masks = block.zero_masks();
        result.logical.screen_points += 64;
        result.logical.screen_batches += 1;
        let packed = masks
            .iter()
            .enumerate()
            .fold(0u64, |v, (g, m)| v | (u64::from(*m) << (16 * g)));
        result.trace = (result.trace ^ step).wrapping_mul(0x100000001b3);
        result.trace = (result.trace ^ packed).wrapping_mul(0x100000001b3);
        if packed == 0 {
            continue;
        }
        result.logical.screen_passes += u64::from(packed.count_ones());
        for r in 0..4 {
            let g = legacy_group64(step, r);
            let mut candidates = masks[g];
            while candidates != 0 {
                let lane = candidates.trailing_zeros() as usize;
                candidates &= candidates - 1;
                let low = g * 16 + lane;
                let full = cursor.value(&offsets, low);
                result.logical.screen_full_checks += 1;
                result.trace = (result.trace ^ low as u64).wrapping_mul(0x100000001b3);
                result.trace = (result.trace ^ u64::from(full)).wrapping_mul(0x100000001b3);
                if full == 0 {
                    result.model = Some(((step ^ (step >> 1)) << 6) | low as u64);
                    result.complete = true;
                    return result;
                }
                result.logical.screen_full_rejected += 1;
            }
        }
    }
    result.complete = true;
    result
}
fn solve_byte(system: &System, n: u8, limit: u64, native: bool) -> Solved {
    if n > 24 {
        return solve_packed(system, n, limit);
    }
    let tick = Instant::now();
    let Some(form) = SyndromeForm::from_system(system, n) else {
        return solve_packed(system, n, limit);
    };
    let cap = if limit == 0 { 0 } else { 1 << 24 };
    let got = if native {
        enumerate_bytes::<NativeByte64>(&form, cap)
    } else {
        enumerate_bytes::<ScalarByte64>(&form, cap)
    };
    Solved {
        outcome: if !got.complete {
            Outcome::Unknown("SCREEN_CAP")
        } else {
            got.model.map_or(Outcome::Unsat, Outcome::Sat)
        },
        logical: got.logical,
        profile: Profile {
            partial_ns: tick.elapsed().as_nanos(),
            kernel_calls_by_active: vec![0; n as usize + 1],
            ..Profile::default()
        },
        trace: got.trace,
    }
}

#[derive(Clone, Copy)]
struct Word64([NativeDeltaBlock; 4]);
#[cfg(target_arch = "aarch64")]
fn native_first_zero16(block: &NativeDeltaBlock) -> Option<u8> {
    unsafe {
        use std::arch::aarch64::*;
        let [a, b, c, d] = block.0;
        if vminvq_u32(vminq_u32(vminq_u32(a, b), vminq_u32(c, d))) != 0 {
            return None;
        }
    }
    block
        .values()
        .iter()
        .flatten()
        .position(|&v| v == 0)
        .map(|v| v as u8)
}
#[cfg(target_arch = "x86_64")]
fn native_first_zero16(block: &NativeDeltaBlock) -> Option<u8> {
    unsafe {
        use std::arch::x86_64::*;
        let [a, b, c, d] = block.0;
        let z = _mm_setzero_si128();
        let any = _mm_or_si128(
            _mm_or_si128(_mm_cmpeq_epi32(a, z), _mm_cmpeq_epi32(b, z)),
            _mm_or_si128(_mm_cmpeq_epi32(c, z), _mm_cmpeq_epi32(d, z)),
        );
        if _mm_movemask_epi8(any) == 0 {
            return None;
        }
    }
    block
        .values()
        .iter()
        .flatten()
        .position(|&v| v == 0)
        .map(|v| v as u8)
}
#[cfg(not(any(target_arch = "aarch64", target_arch = "x86_64")))]
fn native_first_zero16(block: &NativeDeltaBlock) -> Option<u8> {
    block
        .values()
        .iter()
        .flatten()
        .position(|&v| v == 0)
        .map(|v| v as u8)
}
fn add_screen_work(target: &mut Logical, source: &Logical) {
    target.screen_calls += source.screen_calls;
    target.screen_points += source.screen_points;
    target.screen_batches += source.screen_batches;
    target.screen_passes += source.screen_passes;
    target.screen_full_checks += source.screen_full_checks;
    target.screen_full_rejected += source.screen_full_rejected;
    target.screen_fallback_points += source.screen_fallback_points;
    target.screen_fallback_batches += source.screen_fallback_batches;
}
impl Word64 {
    fn from_values(values: &[u32; 64]) -> Self {
        Self(std::array::from_fn(|g| {
            NativeDeltaBlock::from_values(&std::array::from_fn(|r| {
                std::array::from_fn(|j| values[g * 16 + r * 4 + j])
            }))
        }))
    }
    fn affine(constant: u32, linear: &[u32; 6]) -> Self {
        let mut values = [constant; 64];
        for y in 1usize..64 {
            let bit = y & y.wrapping_neg();
            values[y] = values[y ^ bit] ^ linear[bit.trailing_zeros() as usize];
        }
        Self::from_values(&values)
    }
    fn xor_uniform(&mut self, value: u32) {
        for g in &mut self.0 {
            g.xor_uniform(value);
        }
    }
    fn xor_block(&mut self, other: &Self) {
        for i in 0..4 {
            self.0[i].xor_block(&other.0[i]);
        }
    }
}
fn enumerate_word64(form: &SyndromeForm, cap: u64) -> Enumeration {
    if form.n < 6 {
        return enumerate_quiet(form, cap);
    }
    let mut cursor = SixCursor::new(form);
    let offsets = low_quadratic_offsets64(form);
    let mut block = Word64::affine(form.constant, &cursor.linear);
    block.xor_block(&Word64::from_values(&offsets));
    let cross: Vec<_> = cursor.cross.iter().map(|a| Word64::affine(0, a)).collect();
    let mut out = Enumeration {
        model: None,
        complete: false,
        points: 0,
        batches: 0,
        checksum: 0,
    };
    for step in 0..1u64 << (form.n - 6) {
        if cap - out.points < 64 {
            return out;
        }
        if let Some((j, d)) = cursor.advance::<false>(form, step) {
            block.xor_uniform(d);
            block.xor_block(&cross[j]);
        }
        out.points += 64;
        out.batches += 1;
        for r in 0..4 {
            let g = legacy_group64(step, r);
            if let Some(lane) = native_first_zero16(&block.0[g]) {
                out.model = Some(((step ^ (step >> 1)) << 6) | (g as u64 * 16 + u64::from(lane)));
                out.complete = true;
                return out;
            }
        }
    }
    out.complete = true;
    out
}
fn solve_gray_control(system: &System, n: u8, limit: u64, wide: bool) -> Solved {
    if n > 24 {
        return solve_packed(system, n, limit);
    }
    let tick = Instant::now();
    let Some(form) = SyndromeForm::from_system(system, n) else {
        return solve_packed(system, n, limit);
    };
    let cap = if limit == 0 { 0 } else { 1 << 24 };
    let got = if wide {
        enumerate_word64(&form, cap)
    } else {
        enumerate_quiet(&form, cap)
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
        trace: 0,
    }
}
