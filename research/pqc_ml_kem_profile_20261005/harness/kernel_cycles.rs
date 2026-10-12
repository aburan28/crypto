/// Scratch-only kernel timings, for comparison with pq-crystals' test_speed.
#[doc(hidden)]
pub fn kernel_cycles() -> Vec<(&'static str, u64)> {
    fn rd() -> u64 {
        unsafe { core::arch::x86_64::_rdtsc() }
    }
    fn med(mut v: Vec<u64>) -> u64 {
        v.sort_unstable();
        v[v.len() / 2]
    }
    let n = 2000;
    let mut out = vec![];
    let mut f = Poly::zero();
    for (i, c) in f.0.iter_mut().enumerate() {
        *c = (i as i16 * 13) % Q;
    }
    let g = f;
    let time = |name: &'static str, out: &mut Vec<(&'static str, u64)>, mut op: Box<dyn FnMut()>| {
        for _ in 0..100 {
            op();
        }
        let v: Vec<u64> = (0..n)
            .map(|_| {
                let t = rd();
                op();
                rd() - t
            })
            .collect();
        out.push((name, med(v)));
    };
    let mut x = f;
    time("ntt", &mut out, Box::new(move || {
        x = g;
        ntt(&mut x);
        std::hint::black_box(&x);
    }));
    let mut x = f;
    time("invntt", &mut out, Box::new(move || {
        x = g;
        invntt(&mut x);
        std::hint::black_box(&x);
    }));
    let a = [f; KMAX];
    let b = [g; KMAX];
    time("basemul_acc k=3", &mut out, Box::new(move || {
        let mut r = Poly::zero();
        for j in 0..3 {
            poly_basemul_acc(&mut r, &a[j], &b[j]);
        }
        poly_reduce(&mut r);
        std::hint::black_box(&r);
    }));
    let rho = [7u8; 32];
    for (name, k) in [("gen_a k=2", 2usize), ("gen_a k=3", 3), ("gen_a k=4", 4)] {
        time(name, &mut out, Box::new(move || {
            let mut m = [[Poly::zero(); KMAX]; KMAX];
            expand_matrix(&rho, k, &mut m);
            std::hint::black_box(&m);
        }));
    }
    let s = [3u8; 32];
    time("getnoise eta2 x1", &mut out, Box::new(move || {
        std::hint::black_box(sample_cbd_prf(&s, 0, 2));
    }));
    time("getnoise eta2 x7 (batched)", &mut out, Box::new(move || {
        let mut o = [Poly::zero(); 7];
        sample_noise(&s, 0, &[2; 7], &mut o);
        std::hint::black_box(&o);
    }));
    let mut c = [0u8; 320];
    time("compress d=10", &mut out, Box::new(move || {
        poly_compress_encode(&mut c, &g, 10);
        std::hint::black_box(&c);
    }));
    let cb = [0x5au8; 320];
    time("decompress d=10", &mut out, Box::new(move || {
        std::hint::black_box(poly_decode_decompress(&cb, 10));
    }));
    for (name, d) in [("compress d=4", 4usize), ("compress d=1 (tomsg)", 1), ("compress d=11", 11)] {
        let mut c = [0u8; 352];
        time(name, &mut out, Box::new(move || {
            poly_compress_encode(&mut c[..32 * d], &g, d);
            std::hint::black_box(&c);
        }));
    }
    for (name, d) in [("decompress d=4", 4usize), ("decompress d=1 (frommsg)", 1)] {
        let cb = [0x5au8; 352];
        time(name, &mut out, Box::new(move || {
            std::hint::black_box(poly_decode_decompress(&cb[..32 * d], d));
        }));
    }
    let ekb = vec![0x11u8; 384 * 3 + 32];
    time("decode_t_hat k=3 (frombytes)", &mut out, Box::new(move || {
        let mut t = [Poly::zero(); KMAX];
        decode_t_hat(3, &ekb, &mut t);
        std::hint::black_box(&t);
    }));
    let m = [9u8; 32];
    time("sha3_512_2 (G)", &mut out, Box::new(move || {
        std::hint::black_box(sha3_512_2(&m, &[1u8; 32]));
    }));
    time("zero Matrix (KMAX^2 polys)", &mut out, Box::new(move || {
        let a = [[Poly::zero(); KMAX]; KMAX];
        std::hint::black_box(&a);
    }));
    time("vec![0; ct_len 768]", &mut out, Box::new(move || {
        std::hint::black_box(vec![0u8; 1088]);
    }));
    let p = crypto_lib::pqc::ml_kem::ML_KEM_768;
    let (ek768, _) = ml_kem_keygen_internal(&p, &[1u8; 32], &[2u8; 32]);
    let pe = MlKemPreparedEncapsKey::new(&p, &ek768).unwrap();
    let (a768, t768) = (pe.a.clone(), pe.t_hat);
    time("kpke_encrypt_with k=3 (no hashing of ek, no gen_a)", &mut out, Box::new(move || {
        let mut c = [0u8; 1088];
        kpke_encrypt_with(&p, &a768, &t768, &[3u8; 32], &[4u8; 32], &mut c);
        std::hint::black_box(&c);
    }));
    let ek2 = ek768.clone();
    time("encaps_internal k=3 (total)", &mut out, Box::new(move || {
        std::hint::black_box(ml_kem_encaps_internal(&p, &ek2, &[5u8; 32]));
    }));
    let keys: Vec<_> = (0..16u8).map(|i| ml_kem_keygen_internal(&p, &[i; 32], &[2u8; 32]).0).collect();
    let mut idx = 0usize;
    time("encaps_internal k=3 (16 keys rotating)", &mut out, Box::new(move || {
        idx = (idx + 1) % 16;
        std::hint::black_box(ml_kem_encaps_internal(&p, &keys[idx], &[idx as u8; 32]));
    }));
    let ekr = ek768.clone();
    let mut j = 0u8;
    time("encaps_internal k=3 (same key, rotating m)", &mut out, Box::new(move || {
        j = j.wrapping_add(1);
        std::hint::black_box(ml_kem_encaps_internal(&p, &ekr, &[j; 32]));
    }));
    let ek = vec![0x11u8; 384 * 3 + 32];
    time("sha3_256(ek) k=3", &mut out, Box::new(move || {
        std::hint::black_box(sha3_256(&ek));
    }));
    out
}
