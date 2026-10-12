//! Decapsulation on adversarial ciphertexts must never panic, must be
//! deterministic, and must not return the honest key for anything but the
//! honest ciphertext.  Toy16 parameters.
mod common;
use common::testrng::TestRng;
use eigen_seta::cm::{point_from_x, random_point};
use eigen_seta::field::Field;
use eigen_seta::kem::{Ciphertext, EigenSeta};
use eigen_seta::params::{ParamSet, Toy16};
use eigen_seta::weier::WCurve;
use rand_core::RngCore;

type K = EigenSeta<Toy16, 4>;
type F = eigen_seta::fp2::Fp2<<Toy16 as ParamSet<4>>::Cfg, 4>;

fn enc(xs: &[F]) -> Vec<u8> {
    xs.iter().flat_map(|x| x.to_bytes()).collect()
}

/// A well-formed ciphertext (it passes parsing and point lifting) on the
/// Montgomery curve with coefficient `a`: x(S), x(T), x(S - T) for random
/// points S, T.  Returns None if the curve has no usable points.
fn on_curve_ct(a: F, rng: &mut TestRng) -> Option<Vec<u8>> {
    let w = WCurve::from_montgomery(a);
    let s = random_point(&w, rng);
    let t = random_point(&w, rng);
    let st = w.sub(&s, &t);
    if bool::from(s.inf) || bool::from(t.inf) || bool::from(st.inf) {
        return None;
    }
    Some(enc(&[a, s.x, t.x, st.x]))
}

fn rand_f(rng: &mut TestRng) -> F {
    let mut b = vec![0u8; F::byte_len()];
    loop {
        rng.fill_bytes(&mut b);
        b[0] &= 0x7f; // keep under the modulus most of the time
        b[F::byte_len() / 2] &= 0x7f;
        let (x, ok) = F::from_bytes(&b);
        if bool::from(ok) {
            return x;
        }
    }
}

#[test]
fn decaps_never_panics_on_adversarial_ciphertexts() {
    let mut rng = TestRng(0xC0FFEE);
    let (pk, sk) = K::keygen(&mut rng);
    let (ct, k) = K::encaps(&pk, &mut rng).unwrap();
    assert_eq!(K::decaps(&sk, &ct), k);
    let l = ct.bytes.len();
    let fl = F::byte_len();
    let mut cases: Vec<(String, Vec<u8>)> = vec![];
    // wrong lengths and random bytes (parsing rejects most of these)
    cases.push(("empty".into(), vec![]));
    cases.push(("short".into(), ct.bytes[..l - 1].to_vec()));
    let mut long = ct.bytes.clone();
    long.push(0);
    cases.push(("long".into(), long));
    cases.push(("all-ff".into(), vec![0xff; l]));
    cases.push(("all-zero".into(), vec![0; l]));
    for i in 0..8 {
        let mut b = vec![0u8; l];
        rng.fill_bytes(&mut b);
        cases.push((format!("random-bytes-{i}"), b));
    }
    // honest ciphertext with each component flipped or swapped
    for comp in 0..4 {
        let mut b = ct.bytes.clone();
        b[comp * fl + fl - 1] ^= 1;
        cases.push((format!("flip-component-{comp}"), b));
    }
    let mut sw = ct.bytes.clone();
    sw[fl..2 * fl].copy_from_slice(&ct.bytes[2 * fl..3 * fl]);
    sw[2 * fl..3 * fl].copy_from_slice(&ct.bytes[fl..2 * fl]);
    cases.push(("swap-S-T".into(), sw));
    // honest curve, small-order and degenerate points: x = 0 is (0, 0) of order 2
    let a_c = F::from_bytes(&ct.bytes[..fl]).0;
    let wc = WCurve::from_montgomery(a_c);
    let x0 = F::zero();
    if let Some(sp) = point_from_x(&wc, F::from_counter(1)) {
        let p2 = wc.sub(&sp, &eigen_seta::weier::WPoint::new(x0, F::zero()));
        cases.push(("S-order-2".into(), enc(&[a_c, x0, sp.x, p2.x])));
        cases.push(("S-equals-T".into(), enc(&[a_c, sp.x, sp.x, x0])));
    }
    // special and singular curves: A = 0 (j = 1728), A = +-2 (singular), A = 1, 3
    for (name, a) in [
        ("A=0", F::zero()),
        ("A=2", F::from_i64(2)),
        ("A=-2", F::from_i64(-2)),
        ("A=1", F::one()),
        ("A=3", F::from_i64(3)),
    ] {
        if let Some(b) = on_curve_ct(a, &mut rng) {
            cases.push((format!("curve-{name}"), b));
        }
    }
    // random curves (almost all ordinary) with valid points on them
    let mut made = 0;
    while made < 24 {
        let a = rand_f(&mut rng);
        if let Some(b) = on_curve_ct(a, &mut rng) {
            cases.push((format!("random-curve-{made}"), b));
            made += 1;
        }
    }
    // a foreign but honest ciphertext (another key)
    let (pk2, _) = K::keygen(&mut rng);
    cases.push(("foreign".into(), K::encaps(&pk2, &mut rng).unwrap().0.bytes));

    for (name, b) in &cases {
        let c = Ciphertext { bytes: b.clone() };
        let r1 = std::panic::catch_unwind(|| K::decaps(&sk, &c));
        let r1 = r1.unwrap_or_else(|_| panic!("decaps panicked on case {name}"));
        let r2 = K::decaps(&sk, &c);
        assert_eq!(r1, r2, "decaps not deterministic on case {name}");
        assert_ne!(r1, k, "case {name} returned the honest shared secret");
    }
    eprintln!("robustness: {} adversarial cases, no panic", cases.len());
}

/// Encapsulation to an adversarial public key must terminate without
/// panicking: it either refuses (None) or returns some ciphertext.  The
/// canonical-basis search is capped, so no curve can make it spin.
#[test]
fn encaps_terminates_on_adversarial_public_keys() {
    use eigen_seta::kem::PublicKey;
    let mut rng = TestRng(0xBADC0DE);
    let fl = F::byte_len();
    let mut keys: Vec<(String, Vec<u8>)> = vec![
        ("empty".into(), vec![]),
        ("all-zero".into(), vec![0; 4 * fl]),
        ("all-ff".into(), vec![0xff; 4 * fl]),
    ];
    for (name, a) in [
        ("A=0", F::zero()),
        ("A=2", F::from_i64(2)),
        ("A=-2", F::from_i64(-2)),
    ] {
        if let Some(b) = on_curve_ct(a, &mut rng) {
            keys.push((format!("curve-{name}"), b));
        }
    }
    let mut made = 0;
    while made < 8 {
        if let Some(b) = on_curve_ct(rand_f(&mut rng), &mut rng) {
            keys.push((format!("random-curve-{made}"), b));
            made += 1;
        }
    }
    for (name, b) in &keys {
        let pk = PublicKey { bytes: b.clone() };
        let r = std::panic::catch_unwind(|| {
            let mut r2 = TestRng(7);
            K::encaps(&pk, &mut r2).map(|(c, _)| c.bytes.len())
        });
        assert!(r.is_ok(), "encaps panicked on public key {name}");
    }
    eprintln!(
        "robustness: {} adversarial public keys, encaps terminated",
        keys.len()
    );
}
