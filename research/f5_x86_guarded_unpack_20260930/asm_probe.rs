#[repr(transparent)]
pub struct Mono(pub u64);

#[no_mangle]
pub unsafe extern "C" fn checked_decode(
    row_ptr: *const u64,
    row_len: usize,
    cols_ptr: *const u64,
    cols_len: usize,
    out: *mut Mono,
) -> usize {
    let row = unsafe { std::slice::from_raw_parts(row_ptr, row_len) };
    let cols = unsafe { std::slice::from_raw_parts(cols_ptr, cols_len) };
    let mut written = 0usize;
    for (w, &word) in row.iter().enumerate() {
        let mut bits = word;
        while bits != 0 {
            let c = w * 64 + bits.trailing_zeros() as usize;
            bits &= bits - 1;
            unsafe { out.add(written).write(Mono(cols[c])) };
            written += 1;
        }
    }
    written
}

#[no_mangle]
pub unsafe extern "C" fn guarded_decode(
    row_ptr: *const u64,
    row_len: usize,
    cols_ptr: *const u64,
    cols_len: usize,
    out: *mut Mono,
) -> usize {
    let row = unsafe { std::slice::from_raw_parts(row_ptr, row_len) };
    let cols = unsafe { std::slice::from_raw_parts(cols_ptr, cols_len) };
    let words = cols_len.div_ceil(64);
    let tail_bits = cols_len % 64;
    let invalid_tail = tail_bits != 0
        && row_len == words
        && row.last().is_some_and(|&last| last >> tail_bits != 0);
    if row_len > words || invalid_tail {
        return unsafe { checked_decode(row_ptr, row_len, cols_ptr, cols_len, out) };
    }
    let mut written = 0usize;
    for (w, &word) in row.iter().enumerate() {
        let mut bits = word;
        while bits != 0 {
            let c = w * 64 + bits.trailing_zeros() as usize;
            bits &= bits - 1;
            unsafe { out.add(written).write(Mono(*cols.get_unchecked(c))) };
            written += 1;
        }
    }
    written
}
