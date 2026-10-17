// Exact bounded census for PROTOCOL.md. No cryptographic or timing claims.
use std::collections::BTreeMap;

fn prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    let mut d = 2;
    while d * d <= n {
        if n % d == 0 {
            return false;
        }
        d += 1;
    }
    true
}

fn pow_mod(mut base: u64, mut exponent: u64, modulus: u64) -> u64 {
    let mut result = 1;
    while exponent > 0 {
        if exponent & 1 == 1 {
            result = result * base % modulus;
        }
        base = base * base % modulus;
        exponent >>= 1;
    }
    result
}

fn count_direct(p: u64, b: u64) -> u64 {
    let mut count = 1; // identity
    for x in 0..p {
        let rhs = ((x * x % p) * x + b) % p;
        for y in 0..p {
            if y * y % p == rhs {
                count += 1;
            }
        }
    }
    count
}

fn main() {
    let mut by_prime: BTreeMap<u64, (u64, u64)> = BTreeMap::new();
    let mut scanned = 0u64;
    let mut witnesses = 0u64;
    println!("p\tb\tN\ttrace\td\th\trational_kernel_size");
    for p in 7..=199 {
        if !prime(p) || p % 3 != 1 {
            continue;
        }
        for b in 1..p {
            if pow_mod(b, (p - 1) / 2, p) != p - 1 {
                continue;
            }
            scanned += 1;
            let n = count_direct(p, b);
            // Independent character count using a signed accumulator.
            let mut character_count = (p + 1) as i64;
            for x in 0..p {
                let rhs = ((x * x % p) * x + b) % p;
                if rhs != 0 {
                    character_count += if pow_mod(rhs, (p - 1) / 2, p) == 1 {
                        1
                    } else {
                        -1
                    };
                }
            }
            assert_eq!(n as i64, character_count, "count mismatch at p={p}, b={b}");
            let trace = p as i64 + 1 - n as i64;
            assert!(trace * trace <= (4 * p) as i64, "Hasse check p={p}, b={b}");
            let entry = by_prime.entry(p).or_default();
            entry.0 += 1;
            if prime(n) {
                assert_ne!(n, 3);
                assert_ne!(n % 3, 0);
                entry.1 += 1;
                witnesses += 1;
                println!("{p}\t{b}\t{n}\t{trace}\t3\t1\t1");
            }
        }
    }
    eprintln!("scanned={scanned} prime_order_witnesses={witnesses}");
    eprintln!("p\tnonsquare_b\tprime_order_b");
    for (p, (all, prime_order)) in by_prime {
        eprintln!("{p}\t{all}\t{prime_order}");
    }
    assert!(witnesses > 0, "no witness in frozen search range");
}
