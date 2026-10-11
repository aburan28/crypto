//! The Eigen-Séta KEM.
//!
//! * `keygen` is variable time (key generation is assumed to run in a
//!   trusted environment; see SIDE_CHANNELS.md).
//! * `encaps` and `decaps` are constant time in all secret data (the sender's
//!   kernel scalar and masks; the receiver's Kani scalars): fixed public
//!   loop counts, branch-free formulas, constant-time selections.  Branches
//!   in these functions depend only on public data (parameters, the public
//!   key, the ciphertext's validity as curve points).
//! * IND-CCA: Fujisaki–Okamoto with implicit rejection on the deterministic
//!   trapdoor function s -> (E_C, [a]phi(P), [b]phi(Q)), (a, b) = G(s, pk).

use crate::chain::{decaps_chain, first_step_glues, keygen_chain, Pair};
use crate::cm::{add_vartime, eigenvector_2k, mul_vartime, point_from_x, point_of_order, torsion_basis, Omega};
use crate::ct::{Choice, ConditionallySelectable, ConstantTimeEq};
use crate::field::Field;
use crate::fp::Fp;
use crate::fp2::Fp2;
use crate::isog3::chain3;
use crate::mont::{ladder, ladder3pt, recover_y, xadd, xtpl_iter, MontCurve, ProjCurve, XPoint};
use crate::ntheory::{cornacchia, is_probable_prime, modpos, random_below, sqrt_mod_2power, two_squares, v2};
use crate::pairing::{dlog_3k, full_order_3k, weil_sq};
use crate::params::ParamSet;
use crate::scalar::{OddModulus, Scalar};
use crate::weier::{WCurve, WIso, WPoint};
use core::marker::PhantomData;
use num_bigint::{BigInt, BigUint};
use num_traits::{One, Zero};
use rand_core::{CryptoRng, RngCore};
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;
use zeroize::{Zeroize, ZeroizeOnDrop, Zeroizing};

pub const SS_BYTES: usize = 32;

#[derive(Clone, Debug)]
pub struct PublicKey {
    pub bytes: Vec<u8>,
}

#[derive(Clone, Debug)]
pub struct Ciphertext {
    pub bytes: Vec<u8>,
}

/// Secret key.  `lam_s`, `lam_t`, `cp2` are mod 2^E; `c3`, `cp3` mod 3^F3.
#[derive(Clone, Zeroize, ZeroizeOnDrop)]
pub struct SecretKey {
    lam_s: Scalar,
    lam_t: Scalar,
    cp2: Scalar,
    c3: Scalar,
    cp3: Scalar,
    z: [u8; 32],
    #[zeroize(skip)]
    pk: Vec<u8>,
    #[zeroize(skip)]
    pkh: [u8; 32],
    /// Bit length of the secret degree f (public by design: f is drawn from
    /// [f_lo, 2 f_lo)).  f itself is not retained: decapsulation never uses it.
    #[zeroize(skip)]
    f_bits: u64,
}

impl SecretKey {
    pub fn secret_degree_bits(&self) -> u64 {
        self.f_bits
    }
}

/// Redacting Debug: never print secret scalars.
impl core::fmt::Debug for SecretKey {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("SecretKey").field("pkh", &self.pkh).finish_non_exhaustive()
    }
}

fn shake(parts: &[&[u8]], out: usize) -> Vec<u8> {
    let mut h = Shake256::default();
    for p in parts {
        h.update(p);
    }
    let mut r = h.finalize_xof();
    let mut v = vec![0u8; out];
    r.read(&mut v);
    v
}

pub struct EigenSeta<P: ParamSet<N>, const N: usize>(PhantomData<P>);

/// Canonical basis of E_A[3^F3] (public, derived from the public key).
#[derive(Clone, Copy, Debug)]
pub struct PkBasis<Fl: Field> {
    pd: WPoint<Fl>,
    qd: WPoint<Fl>,
    pqd: WPoint<Fl>,
}

/// Public data derived from (pk, ct) by `decaps_prepare`.
#[derive(Clone, Copy, Debug)]
pub struct DecapsPrep<Fl: Field> {
    pkf: [Fl; 4],
    a_c: Fl,
    s_pt: WPoint<Fl>,
    t_pt: WPoint<Fl>,
    x1: WPoint<Fl>,
    x2: WPoint<Fl>,
    pk_basis: PkBasis<Fl>,
}

type F<P, const N: usize> = Fp2<<P as ParamSet<N>>::Cfg, N>;

/// Public parameters derived at run time.
struct Ctx<P: ParamSet<N>, const N: usize> {
    d: BigUint,
    d_bits: usize,
    e: usize,
    e_limbs: usize,
    modd: OddModulus,
    _p: PhantomData<P>,
}

impl<P: ParamSet<N>, const N: usize> Ctx<P, N> {
    fn new() -> Self {
        let d = BigUint::from(3u32).pow(P::F3 as u32);
        let d_bits = d.bits() as usize;
        let modd = OddModulus::new(&d);
        Ctx { d, d_bits, e: P::E, e_limbs: P::E.div_ceil(64) + 1, modd, _p: PhantomData }
    }
}

fn fp_from_hex<P: ParamSet<N>, const N: usize>(h: &str) -> Fp<P::Cfg, N> {
    let b = BigUint::parse_bytes(h.as_bytes(), 16).unwrap();
    let mut l = [0u64; N];
    for (i, x) in b.to_u64_digits().iter().enumerate() {
        l[i] = *x;
    }
    Fp::from_canonical_limbs(l)
}

/// Deterministic basis (P, Q) of E[3^F3] on the Montgomery curve y^2 = x^3 +
/// A x^2 + x, computed from public data only (variable time).
fn canonical_basis_3<P: ParamSet<N>, const N: usize>(a: &F<P, N>) -> (WPoint<F<P, N>>, WPoint<F<P, N>>) {
    let w = WCurve::from_montgomery(*a);
    let d = BigUint::from(3u32).pow(P::F3 as u32);
    let cof = crate::cm::p_plus_1::<F<P, N>>() / &d;
    let d1 = &d / 3u32;
    let mut found: Vec<WPoint<F<P, N>>> = vec![];
    let mut k = 1u64;
    loop {
        let x = <F<P, N> as Field>::from_counter(k);
        k += 1;
        let Some(r) = point_from_x(&w, x) else { continue };
        let r = mul_vartime(&w, &r, &cof);
        let r3 = mul_vartime(&w, &r, &d1);
        if bool::from(r3.inf) {
            continue;
        }
        if let Some(p0) = found.first() {
            let p3 = mul_vartime(&w, p0, &d1);
            if bool::from(p3.x.ct_eq(&r3.x)) {
                continue; // [3^(k-1)] R in <[3^(k-1)] P>
            }
            return (*p0, r);
        }
        found.push(r);
    }
}

/// Full point [k]P on the Montgomery curve (A), constant time in k and P:
/// ladder + Okeya–Sakurai y-recovery + one inversion.
fn smul_ct<Fl: Field>(a: &Fl, p: &WPoint<Fl>, k: &Scalar, nbits: usize) -> WPoint<Fl> {
    let pc = ProjCurve::from_affine(*a);
    let (r0, r1) = ladder(&XPoint::from_affine(p.x), k, nbits, &pc);
    let (x, y, z) = recover_y(a, &p.x, &p.y, &r0, &r1);
    let zi = z.invert();
    let r = WPoint { x: x * zi, y: y * zi, inf: Choice::from(0) };
    WPoint::conditional_select(&r, &WPoint::infinity(), z.is_zero())
}

fn fbytes<Fl: Field>(xs: &[Fl]) -> Vec<u8> {
    let mut v = vec![];
    for x in xs {
        v.extend(x.to_bytes());
    }
    v
}

fn parse4<Fl: Field>(b: &[u8]) -> Option<[Fl; 4]> {
    let l = Fl::byte_len();
    if b.len() != 4 * l {
        return None;
    }
    let mut out = [Fl::zero(); 4];
    for i in 0..4 {
        let (v, ok) = Fl::from_bytes(&b[i * l..(i + 1) * l]);
        if !bool::from(ok) {
            return None;
        }
        out[i] = v;
    }
    Some(out)
}

impl<P: ParamSet<N>, const N: usize> EigenSeta<P, N> {
    pub fn name() -> &'static str {
        P::NAME
    }
    pub fn pk_bytes() -> usize {
        4 * <F<P, N> as Field>::byte_len()
    }
    pub fn ct_bytes() -> usize {
        4 * <F<P, N> as Field>::byte_len()
    }

    fn ed_curve() -> WCurve<F<P, N>> {
        WCurve::new(
            F::<P, N>::ZERO,
            Fp2::from_fp(fp_from_hex::<P, N>(P::ED_A4_HEX)),
            Fp2::from_fp(fp_from_hex::<P, N>(P::ED_A6_HEX)),
        )
    }

    /// Masks (a, b), odd, mod 2^E, from the kernel scalar s.
    fn masks(ctx: &Ctx<P, N>, s_bytes: &[u8], pkh: &[u8]) -> (Scalar, Scalar) {
        let nb = ctx.e.div_ceil(8);
        let h = Zeroizing::new(shake(&[b"EigenSeta-G", s_bytes, pkh], 2 * nb));
        let mut a = Scalar::from_le_bytes(&h[..nb], ctx.e_limbs).mod_2k(ctx.e);
        let mut b = Scalar::from_le_bytes(&h[nb..], ctx.e_limbs).mod_2k(ctx.e);
        a.l[0] |= 1;
        b.l[0] |= 1;
        (a, b)
    }

    fn s_to_bytes(ctx: &Ctx<P, N>, s: &Scalar) -> Vec<u8> {
        s.to_le_bytes(ctx.d_bits.div_ceil(8))
    }

    /// Deterministic encryption of the kernel scalar s.  Constant time in s.
    /// Public, variable-time: the canonical basis of E_A[3^F3] and x(P - Q).
    fn pk_basis(a_a: &F<P, N>) -> PkBasis<F<P, N>> {
        let (pd, qd) = canonical_basis_3::<P, N>(a_a);
        let pqd = WCurve::from_montgomery(*a_a).sub(&pd, &qd);
        PkBasis { pd, qd, pqd }
    }

    fn enc_det(ctx: &Ctx<P, N>, pk: &[F<P, N>; 4], s: &Scalar, pkh: &[u8]) -> Vec<u8> {
        Self::enc_det_core(ctx, pk, &Self::pk_basis(&pk[0]), s, pkh)
    }

    /// Deterministic encryption given the public basis: constant time in `s`.
    fn enc_det_core(ctx: &Ctx<P, N>, pk: &[F<P, N>; 4], basis: &PkBasis<F<P, N>>, s: &Scalar, pkh: &[u8]) -> Vec<u8> {
        let a_a = pk[0];
        let pc = ProjCurve::from_affine(a_a);
        let PkBasis { pd, qd, pqd } = *basis;
        let xk = ladder3pt(
            &XPoint::from_affine(pd.x),
            &XPoint::from_affine(qd.x),
            &XPoint::from_affine(pqd.x),
            s,
            ctx.d_bits,
            &pc,
        );
        let mut pts = [XPoint::from_affine(pk[1]), XPoint::from_affine(pk[2]), XPoint::from_affine(pk[3])];
        let cod = chain3(&pc, &xk, P::F3, &mut pts);
        let a_c = cod.affine_a();
        let cc = ProjCurve::from_affine(a_c);
        let (xp, xq, xpq) = (pts[0], pts[1], pts[2]);
        let xppq = xadd(&xp, &xq, &xpq); // x(phiP + phiQ)
        let (a, b) = Self::masks(ctx, &Zeroizing::new(Self::s_to_bytes(ctx, s)), pkh);
        let (xs, _) = ladder(&xp, &a, ctx.e, &cc);
        let (xt, _) = ladder(&xq, &b, ctx.e, &cc);
        let m = b.mul_mod_2k(&a.inv_mod_2k(ctx.e), ctx.e);
        let wpt = ladder3pt(&xp, &xq, &xppq, &m, ctx.e, &cc); // phiP - [b/a] phiQ
        let (xst, _) = ladder(&wpt, &a, ctx.e, &cc);
        fbytes(&[a_c, xs.affine(), xt.affine(), xst.affine()])
    }

    pub fn encaps<R: RngCore + CryptoRng>(pk: &PublicKey, rng: &mut R) -> Option<(Ciphertext, [u8; SS_BYTES])> {
        let ctx = Ctx::<P, N>::new();
        let pkf = parse4::<F<P, N>>(&pk.bytes)?;
        let pkh = shake(&[b"EigenSeta-P", &pk.bytes], 32);
        let mut buf = vec![0u8; ctx.d_bits.div_ceil(8) + 16];
        rng.fill_bytes(&mut buf);
        let s_raw = Scalar::from_le_bytes(&buf, buf.len().div_ceil(8));
        let s = ctx.modd.to_scalar_width(&ctx.modd.reduce(&s_raw), ctx.modd.limbs());
        buf.zeroize();
        let ct = Self::enc_det(&ctx, &pkf, &s, &pkh);
        let k = Zeroizing::new(shake(&[b"EigenSeta-K", &Zeroizing::new(Self::s_to_bytes(&ctx, &s)), &ct], SS_BYTES));
        let mut out = [0u8; SS_BYTES];
        out.copy_from_slice(&k);
        Some((Ciphertext { bytes: ct }, out))
    }

    pub fn decaps(sk: &SecretKey, ct: &Ciphertext) -> [u8; SS_BYTES] {
        match Self::decaps_prepare(sk, ct) {
            Some(prep) => Self::decaps_core(sk, &prep, ct),
            // public malformation: the implicit-rejection key
            None => {
                let mut out = [0u8; SS_BYTES];
                out.copy_from_slice(&shake(&[b"EigenSeta-R", &sk.z, &ct.bytes], SS_BYTES));
                out
            }
        }
    }

    /// Public phase of decapsulation (variable time, depends only on the public
    /// key and the ciphertext): parse, lift S and T with their relative sign,
    /// and compute the canonical 3-torsion bases of E_C and E_A.
    pub fn decaps_prepare(sk: &SecretKey, ct: &Ciphertext) -> Option<DecapsPrep<F<P, N>>> {
        let [a_c, xs, xt, xst] = parse4::<F<P, N>>(&ct.bytes)?;
        let pkf = parse4::<F<P, N>>(&sk.pk)?;
        let wc = WCurve::from_montgomery(a_c);
        let s_pt = point_from_x(&wc, xs)?;
        let mut t_pt = point_from_x(&wc, xt)?;
        if !bool::from(wc.sub(&s_pt, &t_pt).x.ct_eq(&xst)) {
            t_pt = t_pt.neg();
            if !bool::from(wc.sub(&s_pt, &t_pt).x.ct_eq(&xst)) {
                return None;
            }
        }
        let (x1, x2) = canonical_basis_3::<P, N>(&a_c);
        Some(DecapsPrep { pkf, a_c, s_pt, t_pt, x1, x2, pk_basis: Self::pk_basis(&pkf[0]) })
    }

    /// Secret phase of decapsulation: trapdoor inversion, deterministic
    /// re-encryption, comparison and selection.  Intended to be constant
    /// time in the secret key and in s; kept out of line so that
    /// instruction-count tools can measure it on its own (see
    /// examples/ct_probe.rs and SIDE_CHANNELS.md).
    #[inline(never)]
    pub fn decaps_core(sk: &SecretKey, prep: &DecapsPrep<F<P, N>>, ct: &Ciphertext) -> [u8; SS_BYTES] {
        let reject = Zeroizing::new(shake(&[b"EigenSeta-R", &sk.z, &ct.bytes], SS_BYTES));
        let mut out = [0u8; SS_BYTES];
        out.copy_from_slice(&reject);
        let ctx = Ctx::<P, N>::new();
        let s = Self::recover_s_core(sk, prep); // Scalar: wiped on drop
        let ct2 = Self::enc_det_core(&ctx, &prep.pkf, &prep.pk_basis, &s, &sk.pkh);
        let same = crate::ct::bytes_eq(&ct2, &ct.bytes);
        let s_bytes = Zeroizing::new(Self::s_to_bytes(&ctx, &s));
        let k = Zeroizing::new(shake(&[b"EigenSeta-K", &s_bytes, &ct.bytes], SS_BYTES));
        for i in 0..SS_BYTES {
            out[i] = u8::conditional_select(&out[i], &k[i], same);
        }
        out
    }

    /// The trapdoor inversion: recover s from a ciphertext.  Constant time in
    /// the secret key; returns None only on public malformation of `ct`.
    pub fn recover_s(sk: &SecretKey, ct: &Ciphertext) -> Option<Scalar> {
        Some(Self::recover_s_core(sk, &Self::decaps_prepare(sk, ct)?))
    }

    fn recover_s_core(sk: &SecretKey, prep: &DecapsPrep<F<P, N>>) -> Scalar {
        let ctx = Ctx::<P, N>::new();
        let DecapsPrep { pkf, a_c, s_pt, t_pt, x1, x2, ref pk_basis } = *prep;
        let a_a = pkf[0];
        let wc = WCurve::from_montgomery(a_c);
        // secret: Kani kernel generators (c'S, -lamS S), (c'T, -lamT T)
        let e = ctx.e;
        let k1 = (smul_ct(&a_c, &s_pt, &sk.cp2, e), smul_ct(&a_c, &s_pt, &sk.lam_s, e).neg());
        let k2 = (smul_ct(&a_c, &t_pt, &sk.cp2, e), smul_ct(&a_c, &t_pt, &sk.lam_t, e).neg());
        let o = WPoint::infinity();
        let (fc, imgs) = decaps_chain(&wc, &k1, &k2, e, &[(x1, o), (o, x1), (x2, o)]);
        // identify the codomain factors with E_C: Phi(X,0) = (alpha X, -c'X), Phi(0,X) = (c'X, ...)
        let cpx1 = smul_ct(&a_c, &x1, &sk.cp3, ctx.d_bits);
        let cpx2 = smul_ct(&a_c, &x2, &sk.cp3, ctx.d_bits);
        let sig = WIso::find(&fc, &wc);
        let mut alpha1 = WPoint::infinity();
        let mut alpha2 = WPoint::infinity();
        for swap in [false, true] {
            let pick = |pr: &Pair<F<P, N>>, first: bool| if first ^ swap { pr.0 } else { pr.1 };
            for sc in [sig, sig.negate()] {
                let ok_c = sc.eval(&pick(&imgs[0], false)).ct_eq(&cpx1.neg()) & sc.eval(&pick(&imgs[2], false)).ct_eq(&cpx2.neg());
                for sa in [sig, sig.negate()] {
                    let ok_a = sa.eval(&pick(&imgs[1], true)).ct_eq(&cpx1);
                    let ok = ok_c & ok_a;
                    alpha1 = WPoint::conditional_select(&alpha1, &sa.eval(&pick(&imgs[0], true)), ok);
                    alpha2 = WPoint::conditional_select(&alpha2, &sa.eval(&pick(&imgs[2], true)), ok);
                }
            }
        }
        // psi'(X) = alpha(X) - c X
        let psi1 = wc.sub(&alpha1, &smul_ct(&a_c, &x1, &sk.c3, ctx.d_bits));
        let psi2 = wc.sub(&alpha2, &smul_ct(&a_c, &x2, &sk.c3, ctx.d_bits));
        // generator G of psi'(E_C[d]) (full order), the other image H = k G;
        // k from pairings with a reference point R making e(G, R) full order.
        let k3 = P::F3;
        let d = &ctx.d;
        let pcc = ProjCurve::from_affine(a_c);
        let g1_full = !xtpl_iter(&XPoint::from_affine(psi1.x), k3 - 1, &pcc).is_infinity() & !psi1.inf;
        let g = WPoint::conditional_select(&psi2, &psi1, g1_full);
        let h = WPoint::conditional_select(&psi1, &psi2, g1_full);
        let xg = WPoint::conditional_select(&x2, &x1, g1_full);
        let xh = WPoint::conditional_select(&x1, &x2, g1_full);
        let e_g1 = weil_sq(&wc, &g, &x1, d);
        let e_g2 = weil_sq(&wc, &g, &x2, d);
        let r1_full = full_order_3k(&e_g1, k3);
        let r = WPoint::conditional_select(&x2, &x1, r1_full);
        let base = F::<P, N>::conditional_select(&e_g2, &e_g1, r1_full);
        let tgt = weil_sq(&wc, &h, &r, d);
        let kk = dlog_3k(&base, &tgt, k3);
        let kk = ctx.modd.to_scalar_width(&ctx.modd.reduce(&kk), ctx.modd.limbs());
        // kernel of psi' (= ker phi_dual): X_h - k X_g
        let khat = wc.sub(&xh, &smul_ct(&a_c, &xg, &kk, ctx.d_bits));
        // phi_dual and the images of X1, X2 in E_A
        let mut pts = [XPoint::from_affine(x1.x), XPoint::from_affine(x2.x)];
        let cod = chain3(&pcc, &XPoint::from_affine(khat.x), k3, &mut pts);
        let a_cod = cod.affine_a();
        let neg = a_cod.ct_eq(&(-a_a));
        let pa = ProjCurve::from_affine(a_a);
        let mut imx = [F::<P, N>::zero(); 2];
        for i in 0..2 {
            let x = pts[i].affine();
            imx[i] = F::<P, N>::conditional_select(&x, &(-x), neg);
        }
        let full0 = !xtpl_iter(&XPoint::from_affine(imx[0]), k3 - 1, &pa).is_infinity();
        let xkg = F::<P, N>::conditional_select(&imx[1], &imx[0], full0);
        let wa = WCurve::from_montgomery(a_a);
        let (yk, _) = wa.rhs(&xkg).sqrt();
        let kg = WPoint::new(xkg, yk);
        // s = y/x where Kg = x Pd + y Qd:  e(Kg, Qd) = g^(u x), e(Pd, Kg) = g^(u y)
        let (pd, qd) = (pk_basis.pd, pk_basis.qd);
        let z1 = weil_sq(&wa, &kg, &qd, d);
        let z2 = weil_sq(&wa, &pd, &kg, d);
        let s = dlog_3k(&z1, &z2, k3);
        ctx.modd.to_scalar_width(&ctx.modd.reduce(&s), ctx.modd.limbs())
    }

    /// Key generation (variable time).
    pub fn keygen<R: RngCore + CryptoRng>(rng: &mut R) -> (PublicKey, SecretKey) {
        loop {
            if let Some(kp) = Self::keygen_attempt(rng) {
                return kp;
            }
        }
    }

    fn keygen_attempt<R: RngCore + CryptoRng>(rng: &mut R) -> Option<(PublicKey, SecretKey)> {
        let ctx = Ctx::<P, N>::new();
        let e = ctx.e;
        let d = ctx.d.clone();
        let n2 = BigUint::one() << e;
        let ed = Self::ed_curve();
        let om = Omega::find(&ed, rng);
        let t = om.trace;
        let s_d = BigInt::from(sqrt_mod_2power(&((BigUint::one() << (e + 2)) - 7u32), e as u64 + 2));
        let mu = modpos(&((BigInt::from(t) + &s_d) / 2), &n2);
        let nu = modpos(&((BigInt::from(t) - &s_d) / 2), &n2);
        let w = |x: &WPoint<F<P, N>>| om.eval(x);
        let pd_e = eigenvector_2k(&ed, &w, &nu, e, rng); // omega-eigenvalue mu
        let qd_e = eigenvector_2k(&ed, &w, &mu, e, rng); // omega-eigenvalue nu
        let s_dm = modpos(&s_d, &n2);
        let lo = BigUint::parse_bytes(P::F_LO_HEX.as_bytes(), 16).unwrap();
        let pchar = crate::cm::p_plus_1::<F<P, N>>() - 1u32;
        // 1. secret degree f, Kani scalars (c, c') with the canonical schedule
        let (f, c, cp) = loop {
            let f = &lo + random_below(rng, &lo);
            let f = if (&f % 2u32).is_zero() { f + 1u32 } else { f };
            if (&f % 3u32).is_zero() || (&f % 7u32).is_zero() {
                continue;
            }
            let two_f = BigInt::from(f.clone()) * 2 - BigInt::from(n2.clone());
            if modpos(&two_f, &d).is_zero() {
                continue;
            }
            let m = &n2 - &d * &d * &f * &f * 7u32;
            let Some((x, y)) = two_squares(&m, rng, 1 << 14) else { continue };
            let (c0, cp0) = if (&x % 2u32).is_zero() { (x, y) } else { (y, x) };
            let s0 = (&f * &s_dm) % &n2;
            let ds0 = BigInt::from((&d * &s0) % &n2);
            let mut choice = None;
            for sc in [1i32, -1] {
                for scp in [1i32, -1] {
                    let c = BigInt::from(c0.clone()) * sc;
                    let cp = BigInt::from(cp0.clone()) * scp;
                    let lam_s = &ds0 + &c;
                    if v2(&(&cp - &lam_s)).max(v2(&(&cp + &lam_s))) == 2 {
                        choice = Some((c.clone(), cp.clone()));
                    }
                }
            }
            if let Some((c, cp)) = choice {
                break (f, c, cp);
            }
        };
        let g = &n2 - &f;
        let fg = &f * &g;
        // 2. alpha = x0 + x1 omega + (y0 + y1 omega) pi with n(alpha) = f g
        let nform = |u0: &BigInt, u1: &BigInt| -> BigInt { u0 * u0 + BigInt::from(t) * u0 * u1 + u1 * u1 * 2 };
        let ymax = (&fg / &pchar).sqrt() + 1u32;
        let ymax_i = BigInt::from(ymax.clone());
        let mut alpha_c = None;
        for _ in 0..20000 {
            let y0 = BigInt::from(random_below(rng, &(&ymax * 2u32 + 1u32))) - &ymax_i;
            let y1 = BigInt::from(random_below(rng, &(&ymax * 2u32 + 1u32))) - &ymax_i;
            if y0.is_zero() && y1.is_zero() {
                continue;
            }
            let tt = BigInt::from(fg.clone()) - BigInt::from(pchar.clone()) * nform(&y0, &y1);
            if tt <= BigInt::zero() {
                continue;
            }
            let tt = tt.to_biguint().unwrap();
            let r7 = (&tt % 7u32).to_u64_digits().first().copied().unwrap_or(0);
            if ![1u64, 2, 4].contains(&r7) || !is_probable_prime(&tt, rng, 2) {
                continue;
            }
            // x^2 + 7 y^2 = T  ->  x0 + x1 omega = (x - t y) + 2 y omega
            if let Some((xx, yy)) = cornacchia(7, &tt, rng) {
                let (xx, yy) = (BigInt::from(xx), BigInt::from(yy));
                alpha_c = Some((&xx - BigInt::from(t) * &yy, &yy * 2, y0, y1));
                break;
            }
        }
        let (x0, x1, y0, y1) = alpha_c?;
        let frob = |p: &WPoint<F<P, N>>| WPoint { x: p.x.frob(), y: p.y.frob(), inf: p.inf };
        let mulz = |p: &WPoint<F<P, N>>, k: &BigInt| {
            let r = mul_vartime(&ed, p, &k.magnitude().clone());
            if k.sign() == num_bigint::Sign::Minus {
                r.neg()
            } else {
                r
            }
        };
        let alpha = |p: &WPoint<F<P, N>>| {
            let pp = frob(p);
            let a = add_vartime(&ed, &mulz(p, &x0), &mulz(&om.eval(p), &x1));
            let b = add_vartime(&ed, &mulz(&pp, &y0), &mulz(&om.eval(&pp), &y1));
            add_vartime(&ed, &a, &b)
        };
        // 3. Kani: kernel {([f] R, alpha R)}
        let (r1, r2) = torsion_basis(&ed, 2, e, rng);
        let fz = BigInt::from(f.clone());
        let k1 = (mulz(&r1, &fz), alpha(&r1));
        let k2 = (mulz(&r2, &fz), alpha(&r2));
        if !first_step_glues(&ed, &k1, &k2, e) {
            return None;
        }
        let r3 = point_of_order(&ed, 3, P::F3, rng);
        let (y3a, y3b) = torsion_basis(&ed, 3, P::F3, rng);
        let o = WPoint::infinity();
        let payload = [
            (add_vartime(&ed, &pd_e, &r3), o),
            (add_vartime(&ed, &qd_e, &r3), o),
            (r3, o),
            (y3a, o),
            (y3b, o),
        ];
        let (ea, eb, imgs) = keygen_chain(&ed, &k1, &k2, e, &payload);
        // 4. which factor carries eps (degree f)?  e(eps X, eps Y) = e(X, Y)^f
        let base = weil_sq(&ed, &y3a, &y3b, &d);
        let want = base.pow_public(&(&f % &d).to_u64_digits());
        let e1 = weil_sq(&ea, &imgs[3].0, &imgs[4].0, &d);
        let e2 = weil_sq(&eb, &imgs[3].1, &imgs[4].1, &d);
        let (curve, eps_p, eps_q) = if e1 == want {
            (ea, ea.sub(&imgs[0].0, &imgs[2].0), ea.sub(&imgs[1].0, &imgs[2].0))
        } else if e2 == want {
            (eb, eb.sub(&imgs[0].1, &imgs[2].1), eb.sub(&imgs[1].1, &imgs[2].1))
        } else {
            return None;
        };
        // 5. Montgomery model (B = 1, correct twist) of E_A
        let half = BigUint::one() << (e - 1);
        let t1 = mul_vartime(&curve, &eps_p, &half);
        let t2 = mul_vartime(&curve, &eps_q, &half);
        let roots = [t1.x, t2.x, -curve.a2 - t1.x - t2.x];
        let mut model = None;
        for r in roots {
            let ap = r.mul_small(3) + curve.a2;
            let bp = (r.mul_small(3) + curve.a2.double()) * r + curve.a4;
            let (s, ok) = bp.sqrt();
            if !bool::from(ok) || !bool::from(s.is_square()) {
                continue;
            }
            let si = s.invert();
            let a_m = ap * si;
            // twist check: a random x on the model has [p+1] x = O
            let wm = WCurve::from_montgomery(a_m);
            let rp = crate::cm::random_point(&wm, rng);
            if bool::from(mul_vartime(&wm, &rp, &(&pchar + 1u32)).inf) {
                model = Some((r, si, a_m));
                break;
            }
        }
        let (r0, si, a_m) = model?;
        // 6. random diagonal rescaling of the eigenbasis
        let l1 = random_below(rng, &n2) | BigUint::one();
        let l2 = random_below(rng, &n2) | BigUint::one();
        let pp = mul_vartime(&curve, &eps_p, &l1);
        let qq = mul_vartime(&curve, &eps_q, &l2);
        let pq = curve.sub(&pp, &qq);
        let to_m = |x: F<P, N>| (x - r0) * si;
        let pk_f = [a_m, to_m(pp.x), to_m(qq.x), to_m(pq.x)];
        let pk_bytes = fbytes(&pk_f);
        // 7. secret key
        let s0 = (&f * &s_dm) % &n2;
        let ds0 = BigInt::from((&d * &s0) % &n2);
        let lam_s = modpos(&(&ds0 + &c), &n2);
        let lam_t = modpos(&(&c - &ds0), &n2);
        // Exclude scalars k = -1 mod the point order: x-only y-recovery returns
        // O for [k]P when [k+1]P = O (SIDE_CHANNELS.md R3).  Probability about
        // 2^-E and 3^-F3 per key; checked so the case cannot occur at all.
        let n2m1 = &n2 - 1u32;
        let dm1 = &d - 1u32;
        let cp2v = modpos(&cp, &n2);
        if lam_s == n2m1 || lam_t == n2m1 || cp2v == n2m1 || modpos(&c, &d) == dm1 || modpos(&cp, &d) == dm1 {
            return None;
        }
        let mut z = [0u8; 32];
        rng.fill_bytes(&mut z);
        let pkh = shake(&[b"EigenSeta-P", &pk_bytes], 32);
        let mut pkh_a = [0u8; 32];
        pkh_a.copy_from_slice(&pkh);
        let el = ctx.e_limbs;
        let dl = ctx.modd.limbs();
        let sk = SecretKey {
            lam_s: Scalar::from_biguint(&lam_s, el),
            lam_t: Scalar::from_biguint(&lam_t, el),
            cp2: Scalar::from_biguint(&cp2v, el),
            c3: Scalar::from_biguint(&modpos(&c, &d), dl),
            cp3: Scalar::from_biguint(&modpos(&cp, &d), dl),
            z,
            pk: pk_bytes.clone(),
            pkh: pkh_a,
            f_bits: f.bits(),
        };
        let pk = PublicKey { bytes: pk_bytes };
        // 8. self-test (one encapsulation must decapsulate)
        let (ctx_, k) = Self::encaps(&pk, rng)?;
        if Self::decaps(&sk, &ctx_) != k {
            return None;
        }
        Some((pk, sk))
    }
}

pub(crate) fn _unused(_: MontCurve<Fp2<crate::params::Toy16Fp, 4>>) {}
