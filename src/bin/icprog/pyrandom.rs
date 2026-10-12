//! CPython's `random.Random`, for the resamples a frozen analysis drew:
//! MT19937 seeded through `init_by_array` from an integer seed, as
//! `random_seed` in `_randommodule.c` seeds it, `getrandbits`, and
//! `randrange(n)` by `_randbelow_with_getrandbits`.  A native bootstrap
//! with the same seed therefore draws the same samples.

const N: usize = 624;
const M: usize = 397;
const MATRIX_A: u32 = 0x9908_b0df;
const UPPER_MASK: u32 = 0x8000_0000;
const LOWER_MASK: u32 = 0x7fff_ffff;

pub struct Random {
    mt: [u32; N],
    index: usize,
}

impl Random {
    fn init_genrand(s: u32) -> Random {
        let mut mt = [0u32; N];
        mt[0] = s;
        for i in 1..N {
            mt[i] = 1_812_433_253u32
                .wrapping_mul(mt[i - 1] ^ (mt[i - 1] >> 30))
                .wrapping_add(i as u32);
        }
        Random { mt, index: N }
    }

    fn init_by_array(key: &[u32]) -> Random {
        let mut r = Random::init_genrand(19_650_218);
        let mt = &mut r.mt;
        let (mut i, mut j) = (1usize, 0usize);
        for _ in 0..N.max(key.len()) {
            mt[i] = (mt[i] ^ (mt[i - 1] ^ (mt[i - 1] >> 30)).wrapping_mul(1_664_525))
                .wrapping_add(key[j])
                .wrapping_add(j as u32);
            i += 1;
            j += 1;
            if i >= N {
                mt[0] = mt[N - 1];
                i = 1;
            }
            if j >= key.len() {
                j = 0;
            }
        }
        for _ in 0..N - 1 {
            mt[i] = (mt[i] ^ (mt[i - 1] ^ (mt[i - 1] >> 30)).wrapping_mul(1_566_083_941))
                .wrapping_sub(i as u32);
            i += 1;
            if i >= N {
                mt[0] = mt[N - 1];
                i = 1;
            }
        }
        mt[0] = 0x8000_0000;
        r
    }

    /// `random.Random(seed)` for a nonnegative integer seed: its 32-bit
    /// words, least significant first, and one zero word for zero.
    pub fn new(seed: u128) -> Random {
        let mut key = Vec::new();
        let mut s = seed;
        while s != 0 {
            key.push(s as u32);
            s >>= 32;
        }
        if key.is_empty() {
            key.push(0);
        }
        Random::init_by_array(&key)
    }

    fn genrand_u32(&mut self) -> u32 {
        let mag01 = |y: u32| if y & 1 == 1 { MATRIX_A } else { 0 };
        if self.index >= N {
            let mt = &mut self.mt;
            for k in 0..N - M {
                let y = (mt[k] & UPPER_MASK) | (mt[k + 1] & LOWER_MASK);
                mt[k] = mt[k + M] ^ (y >> 1) ^ mag01(y);
            }
            for k in N - M..N - 1 {
                let y = (mt[k] & UPPER_MASK) | (mt[k + 1] & LOWER_MASK);
                mt[k] = mt[k + M - N] ^ (y >> 1) ^ mag01(y);
            }
            let y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
            mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01(y);
            self.index = 0;
        }
        let mut y = self.mt[self.index];
        self.index += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^= y >> 18;
        y
    }

    /// `getrandbits(k)` for `1 <= k <= 32`, the only widths drawn here.
    pub fn getrandbits(&mut self, k: u32) -> u32 {
        assert!(
            (1..=32).contains(&k),
            "getrandbits is ported for 1..=32 bits"
        );
        self.genrand_u32() >> (32 - k)
    }

    /// `randrange(n)`: `getrandbits(n.bit_length())` until it is below `n`.
    pub fn randbelow(&mut self, n: u32) -> u32 {
        assert!(n > 0, "empty range for randrange()");
        let k = 32 - n.leading_zeros();
        loop {
            let r = self.getrandbits(k);
            if r < n {
                return r;
            }
        }
    }

    /// `random()`: 53 bits from two words, as `random_random` builds them.
    #[cfg(test)]
    fn random(&mut self) -> f64 {
        let a = self.genrand_u32() >> 5;
        let b = self.genrand_u32() >> 6;
        (a as f64 * 67_108_864.0 + b as f64) * (1.0 / 9_007_199_254_740_992.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_reference_generator_and_pythons_seeding() {
        // mt19937ar.c's own check: init_by_array({0x123, 0x234, 0x345, 0x456}).
        let mut r = Random::init_by_array(&[0x123, 0x234, 0x345, 0x456]);
        let first: Vec<u32> = (0..5).map(|_| r.genrand_u32()).collect();
        assert_eq!(
            first,
            [
                1_067_595_299,
                955_945_823,
                477_289_528,
                4_107_218_783,
                4_228_976_476
            ]
        );
        // random.seed(0); random.random() and random.seed(42); random.random().
        assert_eq!(Random::new(0).random(), 0.844_421_851_525_048_1);
        assert_eq!(Random::new(42).random(), 0.639_426_798_457_883_7);
    }
}
