//! Independent reproduction of rand 0.8 StdRng's registered n17 query law.
//! PCG seed expansion, ChaCha12, then the same multiply-high rejection sampler.
pub struct QueryLaw {
    key: [u32; 8],
    counter: u64,
    words: [u32; 16],
    at: usize,
}
impl QueryLaw {
    pub fn new(seed: u64) -> Self {
        let mut state = seed ^ 0x44455343454e5400;
        let mut key = [0u32; 8];
        for word in &mut key {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(11634580027462260723);
            *word = ((((state >> 18) ^ state) >> 27) as u32).rotate_right((state >> 59) as u32);
        }
        Self {
            key,
            counter: 0,
            words: [0; 16],
            at: 16,
        }
    }
    fn refill(&mut self) {
        fn quarter(s: &mut [u32; 16], a: usize, b: usize, c: usize, d: usize) {
            s[a] = s[a].wrapping_add(s[b]);
            s[d] = (s[d] ^ s[a]).rotate_left(16);
            s[c] = s[c].wrapping_add(s[d]);
            s[b] = (s[b] ^ s[c]).rotate_left(12);
            s[a] = s[a].wrapping_add(s[b]);
            s[d] = (s[d] ^ s[a]).rotate_left(8);
            s[c] = s[c].wrapping_add(s[d]);
            s[b] = (s[b] ^ s[c]).rotate_left(7);
        }
        let mut initial = [0u32; 16];
        initial[..4].copy_from_slice(&[0x61707865, 0x3320646e, 0x79622d32, 0x6b206574]);
        initial[4..12].copy_from_slice(&self.key);
        initial[12] = self.counter as u32;
        initial[13] = (self.counter >> 32) as u32;
        let mut s = initial;
        for _ in 0..6 {
            for [a, b, c, d] in [
                [0, 4, 8, 12],
                [1, 5, 9, 13],
                [2, 6, 10, 14],
                [3, 7, 11, 15],
                [0, 5, 10, 15],
                [1, 6, 11, 12],
                [2, 7, 8, 13],
                [3, 4, 9, 14],
            ] {
                quarter(&mut s, a, b, c, d);
            }
        }
        self.words = std::array::from_fn(|i| s[i].wrapping_add(initial[i]));
        self.counter += 1;
        self.at = 0;
    }
    fn word(&mut self) -> u32 {
        if self.at == 16 {
            self.refill()
        }
        let word = self.words[self.at];
        self.at += 1;
        word
    }
    pub fn scalar(&mut self) -> u64 {
        let width = 65586u64;
        let zone = (width << width.leading_zeros()).wrapping_sub(1);
        loop {
            let next = self.word() as u64 | ((self.word() as u64) << 32);
            let product = next as u128 * width as u128;
            if product as u64 <= zone {
                return 1 + (product >> 64) as u64;
            }
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use rand::{rngs::StdRng, Rng, SeedableRng};
    #[test]
    fn independent_law_matches_native_sampler() {
        for seed in [0, 1, 2026100102, 2026100302, u64::MAX] {
            let mut own = QueryLaw::new(seed);
            let mut native = StdRng::seed_from_u64(seed ^ 0x44455343454e5400);
            for _ in 0..128 {
                assert_eq!(own.scalar(), native.gen_range(1u64..65587));
            }
        }
        let mut law = QueryLaw::new(2026100102);
        law.scalar();
        law.scalar();
        assert_eq!((law.scalar(), law.scalar()), (32326, 42888));
    }
}
