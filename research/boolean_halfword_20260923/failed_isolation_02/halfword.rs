// Exact projection to the first sixteen equation coordinates. A projected zero
// is only a candidate; every hit is checked on the original polynomial system.
trait HalfBlock: Copy {
    const POINTS: usize;
    fn from_values(values: &[u16; 64]) -> Self;
    fn values(&self) -> [u16; 64];
    fn xor_uniform(&mut self, value: u16);
    fn xor_block(&mut self, other: &Self);
    fn any_zero(&self) -> bool;
}
#[derive(Clone, Copy)]
struct ScalarHalf<const P: usize>([u16; P]);
impl<const P: usize> HalfBlock for ScalarHalf<P> {
    const POINTS: usize = P;
    #[inline(always)]
    fn from_values(values: &[u16; 64]) -> Self {
        assert!(P == 16 || P == 64);
        Self(std::array::from_fn(|i| values[i]))
    }
    #[inline(always)]
    fn values(&self) -> [u16; 64] {
        let mut out = [0; 64];
        out[..P].copy_from_slice(&self.0);
        out
    }
    #[inline(always)]
    fn xor_uniform(&mut self, value: u16) {
        for item in &mut self.0 {
            *item ^= value;
        }
    }
    #[inline(always)]
    fn xor_block(&mut self, other: &Self) {
        for i in 0..P {
            self.0[i] ^= other.0[i];
        }
    }
    #[inline(always)]
    fn any_zero(&self) -> bool {
        self.0.contains(&0)
    }
}

#[cfg(target_arch = "aarch64")]
#[derive(Clone, Copy)]
struct NativeHalf<const G: usize>([std::arch::aarch64::uint16x8_t; G]);
#[cfg(target_arch = "aarch64")]
impl<const G: usize> HalfBlock for NativeHalf<G> {
    const POINTS: usize = G * 8;
    #[inline(always)]
    fn from_values(values: &[u16; 64]) -> Self {
        assert!(G == 2 || G == 8);
        // Every load addresses a complete eight-u16 group in the supplied array.
        Self(std::array::from_fn(|g| unsafe {
            std::arch::aarch64::vld1q_u16(values.as_ptr().add(g * 8))
        }))
    }
    #[inline(always)]
    fn values(&self) -> [u16; 64] {
        let mut out = [0; 64];
        for g in 0..G {
            unsafe {
                std::arch::aarch64::vst1q_u16(out.as_mut_ptr().add(g * 8), self.0[g]);
            }
        }
        out
    }
    #[inline(always)]
    fn xor_uniform(&mut self, value: u16) {
        unsafe {
            use std::arch::aarch64::*;
            let word = vdupq_n_u16(value);
            for x in &mut self.0 {
                *x = veorq_u16(*x, word);
            }
        }
    }
    #[inline(always)]
    fn xor_block(&mut self, other: &Self) {
        for g in 0..G {
            self.0[g] = unsafe { std::arch::aarch64::veorq_u16(self.0[g], other.0[g]) };
        }
    }
    #[inline(always)]
    fn any_zero(&self) -> bool {
        unsafe {
            use std::arch::aarch64::*;
            let smallest = if G == 2 {
                vminq_u16(self.0[0], self.0[1])
            } else {
                vminq_u16(
                    vminq_u16(
                        vminq_u16(self.0[0], self.0[1]),
                        vminq_u16(self.0[2], self.0[3]),
                    ),
                    vminq_u16(
                        vminq_u16(self.0[4], self.0[5]),
                        vminq_u16(self.0[6], self.0[7]),
                    ),
                )
            };
            vminvq_u16(smallest) == 0
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[derive(Clone, Copy)]
struct NativeHalf<const G: usize>([std::arch::x86_64::__m128i; G]);
#[cfg(target_arch = "x86_64")]
impl<const G: usize> HalfBlock for NativeHalf<G> {
    const POINTS: usize = G * 8;
    #[inline(always)]
    fn from_values(values: &[u16; 64]) -> Self {
        assert!(G == 2 || G == 8);
        Self(std::array::from_fn(|g| unsafe {
            std::arch::x86_64::_mm_loadu_si128(values.as_ptr().add(g * 8).cast())
        }))
    }
    #[inline(always)]
    fn values(&self) -> [u16; 64] {
        let mut out = [0; 64];
        for g in 0..G {
            unsafe {
                std::arch::x86_64::_mm_storeu_si128(out.as_mut_ptr().add(g * 8).cast(), self.0[g]);
            }
        }
        out
    }
    #[inline(always)]
    fn xor_uniform(&mut self, value: u16) {
        unsafe {
            use std::arch::x86_64::*;
            let word = _mm_set1_epi16(value as i16);
            for x in &mut self.0 {
                *x = _mm_xor_si128(*x, word);
            }
        }
    }
    #[inline(always)]
    fn xor_block(&mut self, other: &Self) {
        for g in 0..G {
            self.0[g] = unsafe { std::arch::x86_64::_mm_xor_si128(self.0[g], other.0[g]) };
        }
    }
    #[inline(always)]
    fn any_zero(&self) -> bool {
        unsafe {
            use std::arch::x86_64::*;
            let zero = _mm_setzero_si128();
            let pair = |i, j| {
                _mm_or_si128(
                    _mm_cmpeq_epi16(self.0[i], zero),
                    _mm_cmpeq_epi16(self.0[j], zero),
                )
            };
            let matches = if G == 2 {
                pair(0, 1)
            } else {
                _mm_or_si128(
                    _mm_or_si128(pair(0, 1), pair(2, 3)),
                    _mm_or_si128(pair(4, 5), pair(6, 7)),
                )
            };
            _mm_movemask_epi8(matches) != 0
        }
    }
}

#[cfg(not(any(target_arch = "aarch64", target_arch = "x86_64")))]
#[derive(Clone, Copy)]
struct NativeHalf<const G: usize>([u16; 64]);
#[cfg(not(any(target_arch = "aarch64", target_arch = "x86_64")))]
impl<const G: usize> HalfBlock for NativeHalf<G> {
    const POINTS: usize = G * 8;
    fn from_values(values: &[u16; 64]) -> Self {
        assert!(G == 2 || G == 8);
        Self(*values)
    }
    fn values(&self) -> [u16; 64] {
        self.0
    }
    fn xor_uniform(&mut self, value: u16) {
        for x in &mut self.0[..Self::POINTS] {
            *x ^= value;
        }
    }
    fn xor_block(&mut self, other: &Self) {
        for i in 0..Self::POINTS {
            self.0[i] ^= other.0[i];
        }
    }
    fn any_zero(&self) -> bool {
        self.0[..Self::POINTS].contains(&0)
    }
}

struct HalfForm {
    n: usize,
    constant: u16,
    linear: [u16; MAX],
    quadratic: [[u16; MAX]; MAX],
}
impl HalfForm {
    fn from_system(system: &System, n: u8) -> Option<Self> {
        if n > 24 || system.len() > 32 {
            return None;
        }
        let mut out = Self {
            n: n as usize,
            constant: 0,
            linear: [0; MAX],
            quadratic: [[0; MAX]; MAX],
        };
        for (e, poly) in system.iter().enumerate() {
            let coefficient = if e < 16 { 1u16 << e } else { 0 };
            for &monomial in poly {
                if monomial >= 1u64 << n {
                    return None;
                }
                if monomial == 0 {
                    out.constant ^= coefficient;
                    continue;
                }
                let rest = monomial & (monomial - 1);
                if rest != 0 && rest & (rest - 1) != 0 {
                    return None;
                }
                // Validate omitted equations too, while encoding only the retained
                // coordinates. The original-system check still receives every row.
                if coefficient == 0 {
                    continue;
                }
                let i = monomial.trailing_zeros() as usize;
                if rest == 0 {
                    out.linear[i] ^= coefficient;
                } else {
                    let j = rest.trailing_zeros() as usize;
                    out.quadratic[i][j] ^= coefficient;
                    out.quadratic[j][i] ^= coefficient;
                }
            }
        }
        Some(out)
    }
    fn from_full(form: &SyndromeForm) -> Self {
        Self {
            n: form.n,
            constant: form.constant as u16,
            linear: form.linear.map(|x| x as u16),
            quadratic: form.quadratic.map(|row| row.map(|x| x as u16)),
        }
    }
    fn value(&self, point: u64) -> u16 {
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
}
struct HalfSchedule<B: HalfBlock, const LOW: usize> {
    block: B,
    differences: [u16; MAX],
    cross: Vec<B>,
}
impl<B: HalfBlock, const LOW: usize> HalfSchedule<B, LOW> {
    fn new(form: &HalfForm) -> Self {
        assert!((LOW == 4 || LOW == 6) && B::POINTS == 1 << LOW && form.n >= LOW);
        let mut values = [0; 64];
        values[0] = form.constant;
        for y in 1usize..B::POINTS {
            let bit = y & y.wrapping_neg();
            let i = bit.trailing_zeros() as usize;
            let rest = y ^ bit;
            let mut value = values[rest] ^ form.linear[i];
            let mut bits = rest;
            while bits != 0 {
                let j = bits.trailing_zeros() as usize;
                bits &= bits - 1;
                value ^= form.quadratic[i][j];
            }
            values[y] = value;
        }
        let mut differences = [0; MAX];
        let cross = (0..form.n - LOW)
            .map(|j| {
                differences[j] = form.linear[j + LOW]
                    ^ if j == 0 {
                        0
                    } else {
                        form.quadratic[j + LOW][j + LOW - 1]
                    };
                let mut flat = [0; 64];
                for y in 1usize..B::POINTS {
                    let bit = y & y.wrapping_neg();
                    flat[y] =
                        flat[y ^ bit] ^ form.quadratic[bit.trailing_zeros() as usize][j + LOW];
                }
                B::from_values(&flat)
            })
            .collect();
        Self {
            block: B::from_values(&values),
            differences,
            cross,
        }
    }
    #[inline(always)]
    fn advance(&mut self, form: &HalfForm, step: u64) {
        if step == 0 {
            return;
        }
        let j = step.trailing_zeros() as usize;
        self.block.xor_uniform(self.differences[j]);
        self.block.xor_block(&self.cross[j]);
        for i in 0..j {
            self.differences[i] ^= form.quadratic[i + LOW][j + LOW];
        }
    }
    #[inline(always)]
    fn fixed<const J: usize>(&mut self, form: &HalfForm, cross: &B) {
        self.block.xor_uniform(self.differences[J]);
        self.block.xor_block(cross);
        for i in 0..J {
            self.differences[i] ^= form.quadratic[i + LOW][J + LOW];
        }
    }
}
#[derive(Default, Clone, Debug, PartialEq, Eq)]
struct HalfWork {
    block_width: usize,
    points: u64,
    batches: u64,
    projected_hits: u64,
    original_checks: u64,
    rejected_hits: u64,
    verified_hits: u64,
}
impl HalfWork {
    fn json(&self) -> String {
        format!("{{\"lane_bits\":16,\"block_width\":{},\"points\":{},\"batches\":{},\"projected_hits\":{},\"original_checks\":{},\"rejected_hits\":{},\"verified_hits\":{}}}",self.block_width,self.points,self.batches,self.projected_hits,self.original_checks,self.rejected_hits,self.verified_hits)
    }
}
struct HalfEnumeration {
    model: Option<u64>,
    complete: bool,
    work: HalfWork,
}
#[cold]
#[inline(never)]
fn half_original_satisfied(system: &System, point: u64) -> bool {
    satisfies(system, point)
}
#[inline(always)]
fn half_check(system: &System, point: u64, out: &mut HalfEnumeration) -> bool {
    out.work.projected_hits += 1;
    out.work.original_checks += 1;
    if half_original_satisfied(system, point) {
        out.work.verified_hits += 1;
        out.model = Some(point);
        out.complete = true;
        true
    } else {
        out.work.rejected_hits += 1;
        false
    }
}
#[inline(always)]
fn half_process<B: HalfBlock, const LOW: usize, const COUNT: bool>(
    block: &B,
    system: &System,
    out: &mut HalfEnumeration,
    step: u64,
) -> bool {
    if COUNT {
        out.work.points += B::POINTS as u64;
        out.work.batches += 1;
    }
    if !block.any_zero() {
        return false;
    }
    let values = block.values();
    let high = (step ^ (step >> 1)) << LOW;
    if LOW == 6 {
        for r in 0..4 {
            let group = legacy_group64(step, r);
            for lane in 0..16 {
                let low = group * 16 + lane;
                if values[low] == 0 && half_check(system, high | low as u64, out) {
                    if !COUNT {
                        out.work.batches = step + 1;
                        out.work.points = (step + 1) << LOW;
                    }
                    return true;
                }
            }
        }
    } else {
        for (low, &value) in values[..16].iter().enumerate() {
            if value == 0 && half_check(system, high | low as u64, out) {
                if !COUNT {
                    out.work.batches = step + 1;
                    out.work.points = (step + 1) << LOW;
                }
                return true;
            }
        }
    }
    false
}
fn half_tiny(form: &HalfForm, system: &System, cap: u64) -> HalfEnumeration {
    let low = if form.n < 4 { 0 } else { 4 };
    let size = 1usize << low;
    let mut out = HalfEnumeration {
        model: None,
        complete: false,
        work: HalfWork {
            block_width: size,
            ..HalfWork::default()
        },
    };
    for step in 0..1u64 << (form.n - low) {
        if cap - out.work.points < size as u64 {
            return out;
        }
        out.work.points += size as u64;
        out.work.batches += 1;
        let high = if low == 0 {
            step
        } else {
            (step ^ (step >> 1)) << low
        };
        for y in 0..size {
            let point = high | y as u64;
            if form.value(point) == 0 && half_check(system, point, &mut out) {
                return out;
            }
        }
    }
    out.complete = true;
    out
}
fn enumerate_half<B: HalfBlock, const LOW: usize>(
    form: &HalfForm,
    system: &System,
    cap: u64,
) -> HalfEnumeration {
    if cap >= 1u64 << form.n {
        enumerate_half_inner::<B, LOW, false>(form, system, cap)
    } else {
        enumerate_half_inner::<B, LOW, true>(form, system, cap)
    }
}
fn enumerate_half_inner<B: HalfBlock, const LOW: usize, const CHECKED: bool>(
    form: &HalfForm,
    system: &System,
    cap: u64,
) -> HalfEnumeration {
    if form.n < LOW {
        return half_tiny(form, system, cap);
    }
    let mut schedule = HalfSchedule::<B, LOW>::new(form);
    let mut out = HalfEnumeration {
        model: None,
        complete: false,
        work: HalfWork {
            block_width: B::POINTS,
            ..HalfWork::default()
        },
    };
    if form.n - LOW < 4 {
        for step in 0..1u64 << (form.n - LOW) {
            if CHECKED && cap - out.work.points < B::POINTS as u64 {
                return out;
            }
            schedule.advance(form, step);
            if half_process::<B, LOW, CHECKED>(&schedule.block, system, &mut out, step) {
                return out;
            }
        }
    } else {
        let low: [B; 4] = std::array::from_fn(|j| schedule.cross[j]);
        for base in (0..1u64 << (form.n - LOW)).step_by(16) {
            if CHECKED && cap - out.work.points < B::POINTS as u64 {
                return out;
            }
            schedule.advance(form, base);
            if half_process::<B, LOW, CHECKED>(&schedule.block, system, &mut out, base) {
                return out;
            }
            macro_rules! step {
                ($j:literal,$r:literal) => {{
                    if CHECKED && cap - out.work.points < B::POINTS as u64 {
                        return out;
                    }
                    schedule.fixed::<$j>(form, &low[$j]);
                    if half_process::<B, LOW, CHECKED>(&schedule.block, system, &mut out, base + $r)
                    {
                        return out;
                    }
                }};
            }
            step!(0, 1);
            step!(1, 2);
            step!(0, 3);
            step!(2, 4);
            step!(0, 5);
            step!(1, 6);
            step!(0, 7);
            step!(3, 8);
            step!(0, 9);
            step!(1, 10);
            step!(0, 11);
            step!(2, 12);
            step!(0, 13);
            step!(1, 14);
            step!(0, 15);
        }
    }
    if !CHECKED {
        out.work.points = 1u64 << form.n;
        out.work.batches = 1u64 << (form.n - LOW);
    }
    out.complete = true;
    out
}
fn solve_half(system: &System, n: u8, limit: u64, arm: &str) -> (Solved, HalfWork) {
    if n > 24 {
        return (
            Solved {
                outcome: Outcome::Unknown("HALF_DOMAIN"),
                logical: Logical::default(),
                profile: Profile::default(),
                trace: 0,
            },
            HalfWork::default(),
        );
    }
    let tick = Instant::now();
    let Some(form) = HalfForm::from_system(system, n) else {
        return (
            Solved {
                outcome: Outcome::Unknown("HALF_DOMAIN"),
                logical: Logical::default(),
                profile: Profile::default(),
                trace: 0,
            },
            HalfWork::default(),
        );
    };
    let cap = if limit == 0 { 0 } else { 1 << 24 };
    let got = match arm {
        "half16_scalar" => enumerate_half::<ScalarHalf<16>, 4>(&form, system, cap),
        "half64_scalar" => enumerate_half::<ScalarHalf<64>, 6>(&form, system, cap),
        "half16_native" => enumerate_half::<NativeHalf<2>, 4>(&form, system, cap),
        "half64_native" => enumerate_half::<NativeHalf<8>, 6>(&form, system, cap),
        "half16_eor3" => half_eor3::<2, 4>(&form, system, cap),
        "half64_eor3" => half_eor3::<8, 6>(&form, system, cap),
        "half_dispatch" if n <= 20 => enumerate_half::<NativeHalf<2>, 4>(&form, system, cap),
        "half_dispatch" => enumerate_half::<NativeHalf<8>, 6>(&form, system, cap),
        _ => panic!("unknown half-word arm"),
    };
    let outcome = if !got.complete {
        Outcome::Unknown("HALF_CAP")
    } else {
        got.model.map_or(Outcome::Unsat, Outcome::Sat)
    };
    let logical = Logical {
        enumeration_points: got.work.points,
        enumeration_batches: got.work.batches,
        enumeration_leaves: 1,
        ..Logical::default()
    };
    let profile = Profile {
        enumeration_ns: tick.elapsed().as_nanos(),
        kernel_calls_by_active: vec![0; n as usize + 1],
        ..Profile::default()
    };
    (
        Solved {
            outcome,
            logical,
            profile,
            trace: 0,
        },
        got.work,
    )
}
