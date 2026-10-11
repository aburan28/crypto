//! SHA-256 implemented from scratch per FIPS 180-4.
//!
//! # Algorithm outline
//! 1. Pre-process: pad the message to a multiple of 512 bits.
//! 2. Parse into 512-bit blocks.
//! 3. For each block, expand the 16-word message schedule to 64 words,
//!    then run the 64-round compression function against eight 32-bit
//!    working variables (a–h), mixing in round constants K[0..63].
//! 4. Add the compressed block into the running hash value.
//! 5. Concatenate the final eight 32-bit words as big-endian bytes.

/// The eight initial hash values — first 32 bits of the fractional parts
/// of the square roots of the first 8 primes.
const H0: [u32; 8] = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
];

/// The 64 round constants — first 32 bits of the fractional parts of the
/// cube roots of the first 64 primes.
const K: [u32; 64] = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
];

// ── Bitwise operations ───────────────────────────────────────────────────────

#[inline]
fn rotr32(x: u32, n: u32) -> u32 {
    x.rotate_right(n)
}
#[inline]
fn ch(e: u32, f: u32, g: u32) -> u32 {
    (e & f) ^ (!e & g)
}
#[inline]
fn maj(a: u32, b: u32, c: u32) -> u32 {
    (a & b) ^ (a & c) ^ (b & c)
}
#[inline]
fn sigma0(a: u32) -> u32 {
    rotr32(a, 2) ^ rotr32(a, 13) ^ rotr32(a, 22)
}
#[inline]
fn sigma1(e: u32) -> u32 {
    rotr32(e, 6) ^ rotr32(e, 11) ^ rotr32(e, 25)
}
#[inline]
fn gamma0(x: u32) -> u32 {
    rotr32(x, 7) ^ rotr32(x, 18) ^ (x >> 3)
}
#[inline]
fn gamma1(x: u32) -> u32 {
    rotr32(x, 17) ^ rotr32(x, 19) ^ (x >> 10)
}

// ── Padding ──────────────────────────────────────────────────────────────────

/// Pad a message to a multiple of 64 bytes (512 bits):
///   append bit '1', then zeros, then the 64-bit big-endian message length.
fn pad(msg: &[u8]) -> Vec<u8> {
    let bit_len = (msg.len() as u64).wrapping_mul(8);
    let mut out = msg.to_vec();
    out.push(0x80);
    while out.len() % 64 != 56 {
        out.push(0x00);
    }
    out.extend_from_slice(&bit_len.to_be_bytes());
    out
}

// ── Compression ──────────────────────────────────────────────────────────────

fn compress(state: &mut [u32; 8], block: &[u8]) {
    // Build the 64-word message schedule W from the 16 input words.
    let mut w = [0u32; 64];
    for i in 0..16 {
        w[i] = u32::from_be_bytes(block[i * 4..i * 4 + 4].try_into().unwrap());
    }
    for i in 16..64 {
        w[i] = gamma1(w[i - 2])
            .wrapping_add(w[i - 7])
            .wrapping_add(gamma0(w[i - 15]))
            .wrapping_add(w[i - 16]);
    }

    let mut a = state[0];
    let mut b = state[1];
    let mut c = state[2];
    let mut d = state[3];
    let mut e = state[4];
    let mut f = state[5];
    let mut g = state[6];
    let mut h = state[7];

    for i in 0..64 {
        let t1 = h
            .wrapping_add(sigma1(e))
            .wrapping_add(ch(e, f, g))
            .wrapping_add(K[i])
            .wrapping_add(w[i]);
        let t2 = sigma0(a).wrapping_add(maj(a, b, c));
        h = g;
        g = f;
        f = e;
        e = d.wrapping_add(t1);
        d = c;
        c = b;
        b = a;
        a = t1.wrapping_add(t2);
    }

    state[0] = state[0].wrapping_add(a);
    state[1] = state[1].wrapping_add(b);
    state[2] = state[2].wrapping_add(c);
    state[3] = state[3].wrapping_add(d);
    state[4] = state[4].wrapping_add(e);
    state[5] = state[5].wrapping_add(f);
    state[6] = state[6].wrapping_add(g);
    state[7] = state[7].wrapping_add(h);
}

// ── Public API ───────────────────────────────────────────────────────────────

/// Compute SHA-256 of `data`, returning a 32-byte digest.
pub fn sha256(data: &[u8]) -> [u8; 32] {
    let padded = pad(data);
    let mut state = H0;

    for block in padded.chunks_exact(64) {
        compress(&mut state, block);
    }

    let mut out = [0u8; 32];
    for (i, word) in state.iter().enumerate() {
        out[i * 4..i * 4 + 4].copy_from_slice(&word.to_be_bytes());
    }
    out
}

/// Compute SHA-224 (truncated SHA-256 with different IV).
pub fn sha224(data: &[u8]) -> [u8; 28] {
    const H224: [u32; 8] = [
        0xc1059ed8, 0x367cd507, 0x3070dd17, 0xf70e5939, 0xffc00b31, 0x68581511, 0x64f98fa7,
        0xbefa4fa4,
    ];
    let padded = pad(data);
    let mut state = H224;
    for block in padded.chunks_exact(64) {
        compress(&mut state, block);
    }
    let mut out = [0u8; 28];
    for (i, word) in state[..7].iter().enumerate() {
        out[i * 4..i * 4 + 4].copy_from_slice(&word.to_be_bytes());
    }
    out
}

// ── Incremental hashing ──────────────────────────────────────────────────────

/// SHA-256 over input fed in pieces, for files too large to hold in memory.
///
/// Blocks are compressed with the x86 SHA extensions when the CPU has them
/// and with [`compress`] otherwise; either way the digest is that of
/// [`sha256`], the portable reference this is tested against.
#[derive(Clone)]
pub struct Sha256 {
    state: [u32; 8],
    block: [u8; 64],
    filled: usize,
    length: u64,
}

impl Default for Sha256 {
    fn default() -> Self {
        Self::new()
    }
}

impl Sha256 {
    pub fn new() -> Self {
        Self {
            state: H0,
            block: [0; 64],
            filled: 0,
            length: 0,
        }
    }

    pub fn update(&mut self, mut data: &[u8]) {
        self.length = self.length.wrapping_add(data.len() as u64);
        if self.filled > 0 {
            let take = data.len().min(64 - self.filled);
            self.block[self.filled..self.filled + take].copy_from_slice(&data[..take]);
            self.filled += take;
            data = &data[take..];
            if self.filled < 64 {
                return;
            }
            compress_blocks(&mut self.state, &self.block);
            self.filled = 0;
        }
        let whole = data.len() - data.len() % 64;
        compress_blocks(&mut self.state, &data[..whole]);
        self.filled = data.len() - whole;
        self.block[..self.filled].copy_from_slice(&data[whole..]);
    }

    pub fn finalize(mut self) -> [u8; 32] {
        let mut tail = [0u8; 128];
        tail[..self.filled].copy_from_slice(&self.block[..self.filled]);
        tail[self.filled] = 0x80;
        let end = if self.filled < 56 { 64 } else { 128 };
        tail[end - 8..end].copy_from_slice(&self.length.wrapping_mul(8).to_be_bytes());
        compress_blocks(&mut self.state, &tail[..end]);
        let mut out = [0u8; 32];
        for (i, word) in self.state.iter().enumerate() {
            out[i * 4..i * 4 + 4].copy_from_slice(&word.to_be_bytes());
        }
        out
    }
}

/// Compress `blocks`, a whole number of 64-byte blocks, into `state`.
fn compress_blocks(state: &mut [u32; 8], blocks: &[u8]) {
    debug_assert_eq!(blocks.len() % 64, 0);
    if blocks.is_empty() {
        return;
    }
    #[cfg(target_arch = "x86_64")]
    if shani::available() {
        // SAFETY: the required CPU features were detected at runtime.
        unsafe { shani::compress(state, blocks) };
        return;
    }
    for block in blocks.chunks_exact(64) {
        compress(state, block);
    }
}

/// The compression function on the SHA extensions: `sha256rnds2` runs two
/// rounds on the state held as the register pair (ABEF, CDGH), and
/// `sha256msg1`/`sha256msg2` extend the message schedule four words at a
/// time.
#[cfg(target_arch = "x86_64")]
mod shani {
    use core::arch::x86_64::*;

    use super::K;

    pub(super) fn available() -> bool {
        std::arch::is_x86_feature_detected!("sha")
            && std::arch::is_x86_feature_detected!("sse2")
            && std::arch::is_x86_feature_detected!("ssse3")
            && std::arch::is_x86_feature_detected!("sse4.1")
    }

    /// # Safety
    ///
    /// The CPU must support SHA, SSE2, SSSE3 and SSE4.1, and `blocks` must
    /// be a whole number of 64-byte blocks.
    #[target_feature(enable = "sha,sse2,ssse3,sse4.1")]
    pub(super) unsafe fn compress(state: &mut [u32; 8], blocks: &[u8]) {
        // Byte-reverses each 32-bit lane: message words are big-endian.
        let bswap = _mm_set_epi64x(0x0c0d_0e0f_0809_0a0b, 0x0405_0607_0001_0203);

        let dcba = _mm_loadu_si128(state.as_ptr().cast());
        let hgfe = _mm_loadu_si128(state.as_ptr().add(4).cast());
        let cdab = _mm_shuffle_epi32(dcba, 0xb1);
        let efgh = _mm_shuffle_epi32(hgfe, 0x1b);
        let mut abef = _mm_alignr_epi8(cdab, efgh, 8);
        let mut cdgh = _mm_blend_epi16(efgh, cdab, 0xf0);

        // Rounds 4i..4i+3 on the schedule words W[4i..4i+4] held in `$w`.
        macro_rules! rounds4 {
            ($w:expr, $i:expr) => {{
                let k = _mm_loadu_si128(K.as_ptr().add(4 * $i).cast());
                let wk = _mm_add_epi32($w, k);
                cdgh = _mm_sha256rnds2_epu32(cdgh, abef, wk);
                abef = _mm_sha256rnds2_epu32(abef, cdgh, _mm_shuffle_epi32(wk, 0x0e));
            }};
        }
        // W[4i..4i+4] from the sixteen words before them, written over the
        // oldest four, then its rounds.
        macro_rules! schedule_rounds4 {
            ($w0:ident, $w1:ident, $w2:ident, $w3:ident, $i:expr) => {{
                let w7 = _mm_alignr_epi8($w3, $w2, 4);
                $w0 = _mm_sha256msg2_epu32(_mm_add_epi32(_mm_sha256msg1_epu32($w0, $w1), w7), $w3);
                rounds4!($w0, $i);
            }};
        }

        for block in blocks.chunks_exact(64) {
            let (abef_in, cdgh_in) = (abef, cdgh);
            let p = block.as_ptr();
            let mut w0 = _mm_shuffle_epi8(_mm_loadu_si128(p.cast()), bswap);
            let mut w1 = _mm_shuffle_epi8(_mm_loadu_si128(p.add(16).cast()), bswap);
            let mut w2 = _mm_shuffle_epi8(_mm_loadu_si128(p.add(32).cast()), bswap);
            let mut w3 = _mm_shuffle_epi8(_mm_loadu_si128(p.add(48).cast()), bswap);

            rounds4!(w0, 0);
            rounds4!(w1, 1);
            rounds4!(w2, 2);
            rounds4!(w3, 3);
            schedule_rounds4!(w0, w1, w2, w3, 4);
            schedule_rounds4!(w1, w2, w3, w0, 5);
            schedule_rounds4!(w2, w3, w0, w1, 6);
            schedule_rounds4!(w3, w0, w1, w2, 7);
            schedule_rounds4!(w0, w1, w2, w3, 8);
            schedule_rounds4!(w1, w2, w3, w0, 9);
            schedule_rounds4!(w2, w3, w0, w1, 10);
            schedule_rounds4!(w3, w0, w1, w2, 11);
            schedule_rounds4!(w0, w1, w2, w3, 12);
            schedule_rounds4!(w1, w2, w3, w0, 13);
            schedule_rounds4!(w2, w3, w0, w1, 14);
            schedule_rounds4!(w3, w0, w1, w2, 15);

            abef = _mm_add_epi32(abef, abef_in);
            cdgh = _mm_add_epi32(cdgh, cdgh_in);
        }

        let feba = _mm_shuffle_epi32(abef, 0x1b);
        let dchg = _mm_shuffle_epi32(cdgh, 0xb1);
        _mm_storeu_si128(state.as_mut_ptr().cast(), _mm_blend_epi16(feba, dchg, 0xf0));
        _mm_storeu_si128(
            state.as_mut_ptr().add(4).cast(),
            _mm_alignr_epi8(dchg, feba, 8),
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn h(hex_str: &str) -> Vec<u8> {
        hex::decode(hex_str).unwrap()
    }

    // ── SHA-256 known-answer tests (FIPS 180-4 + RFC 6234) ────────────────────

    #[test]
    fn sha256_empty() {
        assert_eq!(
            sha256(b"").as_slice(),
            h("e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855").as_slice(),
        );
    }

    #[test]
    fn sha256_abc() {
        // FIPS 180-4 §B.1
        assert_eq!(
            sha256(b"abc").as_slice(),
            h("ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad").as_slice(),
        );
    }

    #[test]
    fn sha256_hello_world() {
        assert_eq!(
            sha256(b"hello world").as_slice(),
            h("b94d27b9934d3e08a52e52d7da7dabfac484efe37a5380ee9088f7ace2efcde9").as_slice(),
        );
    }

    #[test]
    fn sha256_fips180_56byte() {
        // FIPS 180-4 §B.2
        let msg = b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq";
        assert_eq!(
            sha256(msg).as_slice(),
            h("248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1").as_slice(),
        );
    }

    #[test]
    fn sha256_fips180_112byte() {
        // FIPS 180-4 §B.3 — spans 3 blocks
        let msg = b"abcdefghbcdefghicdefghijdefghijkefghijklfghijklmghijklmnhijklmnoijklmnopjklmnopqklmnopqrlmnopqrsmnopqrstnopqrstu";
        assert_eq!(
            sha256(msg).as_slice(),
            h("cf5b16a778af8380036ce59e7b0492370b249b11e8f07a51afac45037afee9d1").as_slice(),
        );
    }

    #[test]
    fn sha256_million_a() {
        // FIPS 180-4 §B.3 — 1,000,000 'a' characters
        let msg = vec![b'a'; 1_000_000];
        assert_eq!(
            sha256(&msg).as_slice(),
            h("cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0").as_slice(),
        );
    }

    #[test]
    fn sha256_block_boundary_55() {
        // 55 bytes — fits exactly in one padded block (1 byte before length).
        // Verified with `printf 'A...A' | shasum -a 256`.
        assert_eq!(
            sha256(&[b'A'; 55]).as_slice(),
            h("8963cc0afd622cc7574ac2011f93a3059b3d65548a77542a1559e3d202e6ab00").as_slice(),
        );
    }

    #[test]
    fn sha256_two_blocks_64() {
        // 64 bytes — forces a second padding block.
        assert_eq!(
            sha256(&[b'B'; 64]).as_slice(),
            h("c422e7070cb1cb455b5de9afee0d975e303d0239c72030cd7414ab5c382d3ae8").as_slice(),
        );
    }

    // ── SHA-224 known-answer tests (FIPS 180-4) ───────────────────────────────

    #[test]
    fn sha224_empty() {
        assert_eq!(
            sha224(b"").as_slice(),
            h("d14a028c2a3a2bc9476102bb288234c415a2b01f828ea62ac5b3e42f").as_slice(),
        );
    }

    #[test]
    fn sha224_abc() {
        // FIPS 180-4 §B.1 (SHA-224)
        assert_eq!(
            sha224(b"abc").as_slice(),
            h("23097d223405d8228642a477bda255b32aadbce4bda0b3f7e36c9da7").as_slice(),
        );
    }

    // ── Incremental hashing ───────────────────────────────────────────────────

    fn next(x: &mut u64) -> u64 {
        *x ^= *x << 13;
        *x ^= *x >> 7;
        *x ^= *x << 17;
        *x
    }

    fn noise(x: &mut u64, len: usize) -> Vec<u8> {
        (0..len).map(|_| next(x) as u8).collect()
    }

    #[test]
    fn streaming_matches_one_shot_at_every_split() {
        let mut x = 0x9E37_79B9_7F4A_7C15u64;
        for len in (0..=200).chain([255, 256, 257, 4095, 4096, 4097, 65_599]) {
            let data = noise(&mut x, len);
            let want = sha256(&data);

            let mut whole = Sha256::new();
            whole.update(&data);
            assert_eq!(whole.finalize(), want, "{len} bytes in one update");

            let mut bytewise = Sha256::new();
            for byte in &data {
                bytewise.update(std::slice::from_ref(byte));
            }
            assert_eq!(bytewise.finalize(), want, "{len} bytes one at a time");

            let mut ragged = Sha256::new();
            let mut rest = &data[..];
            while !rest.is_empty() {
                let n = (1 + next(&mut x) as usize % 150).min(rest.len());
                ragged.update(&rest[..n]);
                ragged.update(&[]);
                rest = &rest[n..];
            }
            assert_eq!(ragged.finalize(), want, "{len} bytes in ragged pieces");
        }
    }

    #[test]
    fn streaming_million_a() {
        let mut hasher = Sha256::default();
        for _ in 0..1000 {
            hasher.update(&[b'a'; 1000]);
        }
        assert_eq!(
            hasher.finalize().as_slice(),
            h("cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0").as_slice(),
        );
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn sha_extensions_agree_with_the_portable_rounds() {
        if !shani::available() {
            eprintln!("no SHA extensions on this CPU: the portable rounds are the only path");
            return;
        }
        let mut x = 0x2545_F491_4F6C_DD1Du64;
        for n in [1usize, 2, 3, 17] {
            let blocks = noise(&mut x, 64 * n);
            let mut state: [u32; 8] = std::array::from_fn(|_| next(&mut x) as u32);
            let mut portable = state;
            for block in blocks.chunks_exact(64) {
                compress(&mut portable, block);
            }
            // SAFETY: the CPU features were detected above.
            unsafe { shani::compress(&mut state, &blocks) };
            assert_eq!(state, portable, "{n} blocks");
        }
    }
}
