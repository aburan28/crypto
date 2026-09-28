"""Source transformations for the next round; no measurements or winner claims."""
import re


LIFT = '''
/// Exact odd-degree point lifting for factor-base construction. General callers,
/// target generation and rho retain the independent general implementation.
fn factor_base_points_with_x(curve: &BinaryCurve, x: &F2mElement) -> Vec<BinaryPoint> {
    if curve.m % 2 == 0 {
        return points_with_x(curve, x);
    }
    let Some(fc) = FastCurve::new(curve) else {
        return points_with_x(curve, x);
    };
    let f = &fc.field;
    let xw = f.from_element(x);
    if xw == 0 {
        return vec![fc.lower(FastPoint::affine(0, f.sqr_k(fc.b, curve.m - 1)))];
    }
    let inverse = f.inv(xw);
    let rhs = xw ^ fc.a ^ f.mul(fc.b, f.sqr(inverse));
    let mut acc = rhs;
    let mut u = rhs;
    for _ in 0..(curve.m - 1) / 2 {
        acc = f.sqr(f.sqr(acc));
        u ^= acc;
    }
    // Checking the equation also rejects trace-one right-hand sides.
    if f.sqr(u) ^ u != rhs { return Vec::new(); }
    let p = FastPoint::affine(xw, f.mul(xw, u));
    vec![fc.lower(p), fc.lower(fc.neg(p))]
}

#[cfg(test)]
mod tournament_lift_equivalence {
    use super::*;
    #[test]
    fn fast_factor_base_lift_matches_general_including_order() {
        let mut checked = 0;
        for (n, a) in [(5,0), (7,1), (13,0), (17,1), (19,0), (19,1), (23,0), (31,1)] {
            let Some(kc) = KoblitzCurve::new(a, n) else { continue; };
            let mask = (1u64 << n) - 1;
            let count = if n <= 13 { 1u64 << n } else { 512 };
            for i in 0..count {
                let raw = if n <= 13 { i } else { i.wrapping_mul(0x9e3779b97f4a7c15) & mask };
                let x = F2mElement::from_biguint(&BigUint::from(raw), n);
                assert_eq!(factor_base_points_with_x(&kc.curve, &x), points_with_x(&kc.curve, &x),
                           "degree {n}, a={a}, x={raw}");
                checked += 1;
            }
        }
        assert!(checked > 8000);
    }
}
'''


def fast_lift(text):
    at = text.index('/// The scalar `λ`')
    text = text[:at] + LIFT + '\n' + text[at:]
    before = 'let Some(point) = points_with_x(&kc.curve, &x).into_iter().next()'
    assert text.count(before) == 1
    text = text.replace(before, 'let Some(point) = factor_base_points_with_x(&kc.curve, &x).into_iter().next()')
    a = text.index('fn finish_factor_base_domain(')
    i = text.index('for p in points_with_x(&kc.curve, x)', a)
    text = text[:i] + text[i:].replace('for p in points_with_x(&kc.curve, x)',
                                     'for p in factor_base_points_with_x(&kc.curve, x)', 1)
    return text


def folded(text):
    a = text.index('    pub fn build(kc: &KoblitzCurve, fb: &FrobeniusFactorBase)')
    i = text.index('Self::build_within(kc, fb, Self::DEFAULT_BYTE_BUDGET)', a)
    text = text[:i] + text[i:].replace('Self::build_within(kc, fb, Self::DEFAULT_BYTE_BUDGET)',
                                     'Self::build_folded_within(kc, fb, Self::DEFAULT_BYTE_BUDGET)', 1)
    a = text.index('    fn build_folded_within(')
    b = a + 1 + re.search(r'\n    (?:pub )?fn ', text[a+1:]).start()
    block = text[a:b]
    assert block.count('.into_par_iter()') == 2
    text = text[:a] + block.replace('.into_par_iter()', '.into_iter()') + text[b:]
    return text
