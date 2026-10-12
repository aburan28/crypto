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
    let ek = vec![0x11u8; 384 * 3 + 32];
    time("sha3_256(ek) k=3", &mut out, Box::new(move || {
        std::hint::black_box(sha3_256(&ek));
    }));
    out
}
