// Same eight-equation projection, transposed into eight 64-point bit planes.
const fn variable_planes() -> [u64; 6] {
    let mut out = [0; 6];
    let mut i = 0;
    while i < 6 {
        let mut y = 0;
        while y < 64 {
            if y & (1 << i) != 0 {
                out[i] |= 1u64 << y;
            }
            y += 1;
        }
        i += 1;
    }
    out
}
const VARIABLE_PLANES: [u64; 6] = variable_planes();
fn planes_from_values(values: &[u8; 64]) -> [u64; 8] {
    let mut out = [0; 8];
    for (y, &v) in values.iter().enumerate() {
        for e in 0..8 {
            out[e] |= u64::from((v >> e) & 1) << y;
        }
    }
    out
}
fn values_from_planes(planes: &[u64; 8]) -> [u8; 64] {
    std::array::from_fn(|y| {
        planes
            .iter()
            .enumerate()
            .fold(0u8, |v, (e, p)| v | (((p >> y) & 1) as u8) << e)
    })
}
#[cfg(target_arch = "aarch64")]
#[derive(Clone, Copy)]
struct NativePlane64([std::arch::aarch64::uint64x2_t; 4]);
#[cfg(target_arch = "aarch64")]
impl ByteBlock64 for NativePlane64 {
    fn from_values(values: &[u8; 64]) -> Self {
        unsafe {
            let planes = planes_from_values(values);
            Self(std::array::from_fn(|g| {
                std::arch::aarch64::vld1q_u64(planes[g * 2..].as_ptr())
            }))
        }
    }
    fn affine(constant: u8, linear: &[u8; 6]) -> Self {
        unsafe {
            use std::arch::aarch64::*;
            Self(std::array::from_fn(|g| {
                let weights = [1u64 << (2 * g), 1u64 << (2 * g + 1)];
                let weights = vld1q_u64(weights.as_ptr());
                let mut value = vtstq_u64(vdupq_n_u64(u64::from(constant)), weights);
                for i in 0..6 {
                    value = veorq_u64(
                        value,
                        vandq_u64(
                            vtstq_u64(vdupq_n_u64(u64::from(linear[i])), weights),
                            vdupq_n_u64(VARIABLE_PLANES[i]),
                        ),
                    );
                }
                value
            }))
        }
    }
    fn values(&self) -> [u8; 64] {
        let mut planes = [0; 8];
        for g in 0..4 {
            unsafe {
                std::arch::aarch64::vst1q_u64(planes[g * 2..].as_mut_ptr(), self.0[g]);
            }
        }
        values_from_planes(&planes)
    }
    fn xor_uniform(&mut self, value: u8) {
        unsafe {
            use std::arch::aarch64::*;
            let v = vdupq_n_u64(u64::from(value));
            for g in 0..4 {
                let weights = [1u64 << (2 * g), 1u64 << (2 * g + 1)];
                self.0[g] = veorq_u64(self.0[g], vtstq_u64(v, vld1q_u64(weights.as_ptr())));
            }
        }
    }
    fn xor_block(&mut self, other: &Self) {
        for g in 0..4 {
            self.0[g] = unsafe { std::arch::aarch64::veorq_u64(self.0[g], other.0[g]) };
        }
    }
    fn zero_masks(&self) -> [u16; 4] {
        unsafe {
            use std::arch::aarch64::*;
            let p = vorrq_u64(
                vorrq_u64(self.0[0], self.0[1]),
                vorrq_u64(self.0[2], self.0[3]),
            );
            let survivors = !(vgetq_lane_u64::<0>(p) | vgetq_lane_u64::<1>(p));
            std::array::from_fn(|g| (survivors >> (16 * g)) as u16)
        }
    }
}
#[cfg(not(target_arch = "aarch64"))]
#[derive(Clone, Copy)]
struct NativePlane64([u64; 8]);
#[cfg(not(target_arch = "aarch64"))]
impl ByteBlock64 for NativePlane64 {
    fn from_values(values: &[u8; 64]) -> Self {
        Self(planes_from_values(values))
    }
    fn affine(constant: u8, linear: &[u8; 6]) -> Self {
        Self(std::array::from_fn(|e| {
            linear.iter().enumerate().fold(
                0u64.wrapping_sub(u64::from((constant >> e) & 1)),
                |v, (i, &a)| v ^ (VARIABLE_PLANES[i] & 0u64.wrapping_sub(u64::from((a >> e) & 1))),
            )
        }))
    }
    fn values(&self) -> [u8; 64] {
        values_from_planes(&self.0)
    }
    fn xor_uniform(&mut self, value: u8) {
        for e in 0..8 {
            self.0[e] ^= 0u64.wrapping_sub(u64::from((value >> e) & 1));
        }
    }
    fn xor_block(&mut self, other: &Self) {
        for e in 0..8 {
            self.0[e] ^= other.0[e];
        }
    }
    fn zero_masks(&self) -> [u16; 4] {
        let survivors = !self.0.iter().fold(0u64, |v, p| v | p);
        std::array::from_fn(|g| (survivors >> (16 * g)) as u16)
    }
}
fn solve_plane(system: &System, n: u8, limit: u64) -> Solved {
    if n > 24 {
        return solve_packed(system, n, limit);
    }
    let tick = Instant::now();
    let Some(form) = SyndromeForm::from_system(system, n) else {
        return solve_packed(system, n, limit);
    };
    let got = enumerate_tiered::<NativePlane64, false>(&form, if limit == 0 { 0 } else { 1 << 24 });
    Solved {
        outcome: if !got.complete {
            Outcome::Unknown("SCREEN_CAP")
        } else {
            got.model.map_or(Outcome::Unsat, Outcome::Sat)
        },
        logical: got.logical,
        profile: Profile {
            partial_ns: tick.elapsed().as_nanos(),
            full_update_words: got.linear_updates,
            secondary_update_words: got.secondary_updates,
            kernel_calls_by_active: vec![0; n as usize + 1],
            ..Profile::default()
        },
        trace: 0,
    }
}
