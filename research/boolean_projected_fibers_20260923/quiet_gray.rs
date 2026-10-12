// Accounting control: identical 16-point ordering without enumeration checksums.
#[cfg(target_arch = "aarch64")]
fn quiet_syndrome_block(constant: u32, linear: &[u32; 4], offsets: &[[u32; 4]; 4]) -> Option<u8> {
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

        let any = vminvq_u32(vminq_u32(vminq_u32(a, b), vminq_u32(c, d))) == 0;
        if any {
            scalar_syndrome_block(constant, linear, offsets).1
        } else {
            None
        }
    }
}
#[cfg(target_arch = "x86_64")]
fn quiet_syndrome_block(constant: u32, linear: &[u32; 4], offsets: &[[u32; 4]; 4]) -> Option<u8> {
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
        let zero = _mm_setzero_si128();
        let matches = _mm_or_si128(
            _mm_or_si128(_mm_cmpeq_epi32(a, zero), _mm_cmpeq_epi32(b, zero)),
            _mm_or_si128(_mm_cmpeq_epi32(c, zero), _mm_cmpeq_epi32(d, zero)),
        );
        if _mm_movemask_epi8(matches) != 0 {
            scalar_syndrome_block(constant, linear, offsets).1
        } else {
            None
        }
    }
}
#[cfg(not(any(target_arch = "aarch64", target_arch = "x86_64")))]
fn quiet_syndrome_block(constant: u32, linear: &[u32; 4], offsets: &[[u32; 4]; 4]) -> Option<u8> {
    scalar_syndrome_block(constant, linear, offsets).1
}

fn enumerate_quiet(form: &SyndromeForm, point_cap: u64) -> Enumeration {
    let mut result = Enumeration {
        model: None,
        complete: false,
        points: 0,
        batches: 0,
        checksum: 0,
    };
    if form.n < 4 {
        for point in 0..1u64 << form.n {
            if result.points == point_cap {
                return result;
            }
            let value = form.value(point);
            result.points += 1;
            result.batches += 1;
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
        let found = quiet_syndrome_block(cursor.constant, &cursor.linear, &offsets);
        result.points += 16;
        result.batches += 1;
        if let Some(low) = found {
            result.model = Some(((step ^ (step >> 1)) << 4) | u64::from(low));
            result.complete = true;
            return result;
        }
    }
    result.complete = true;
    result
}
