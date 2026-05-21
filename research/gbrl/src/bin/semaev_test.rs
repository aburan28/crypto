// Validate semaev_s3: pick concrete points P1, P2 on a short-Weierstrass curve
// over F_65521, compute P3 = -(P1+P2), and verify that S_3(x1, x2, x3) evaluates
// to zero. Also check (a) symmetry under coordinate permutation, (b) that S_3 is
// non-zero on a triple that does NOT sum to O.
//
//     cargo run --release --bin semaev_test

use gbrl::field::{Fp, P};
use gbrl::semaev::{semaev_s2, semaev_s3, semaev_s4, Curve};

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Point {
    Infinity,
    Affine(Fp, Fp),
}

fn on_curve(curve: &Curve, x: Fp, y: Fp) -> bool {
    y * y == x * x * x + curve.a * x + curve.b
}

fn neg(p: Point) -> Point {
    match p {
        Point::Infinity => Point::Infinity,
        Point::Affine(x, y) => Point::Affine(x, -y),
    }
}

fn add(curve: &Curve, p: Point, q: Point) -> Point {
    match (p, q) {
        (Point::Infinity, _) => q,
        (_, Point::Infinity) => p,
        (Point::Affine(x1, y1), Point::Affine(x2, y2)) => {
            if x1 == x2 {
                if y1 + y2 == Fp::zero() {
                    return Point::Infinity;
                }
                let two = Fp::new(2);
                let three = Fp::new(3);
                let m = (three * x1 * x1 + curve.a) * (two * y1).inv().unwrap();
                let x3 = m * m - x1 - x1;
                let y3 = m * (x1 - x3) - y1;
                Point::Affine(x3, y3)
            } else {
                let m = (y2 - y1) * (x2 - x1).inv().unwrap();
                let x3 = m * m - x1 - x2;
                let y3 = m * (x1 - x3) - y1;
                Point::Affine(x3, y3)
            }
        }
    }
}

// Brute-force search for an affine point with x in a small range. Fine for
// validation; swap for Tonelli-Shanks once we move to bigger fields.
fn find_point(curve: &Curve, start_x: u64) -> Point {
    for x_i in start_x..start_x + 1000 {
        let x = Fp(x_i % P);
        let rhs = x * x * x + curve.a * x + curve.b;
        // y = sqrt(rhs) via brute force; fine for tiny P.
        for y_i in 0..P {
            let y = Fp(y_i);
            if y * y == rhs {
                return Point::Affine(x, y);
            }
        }
    }
    panic!("no affine point found in x range starting at {}", start_x);
}

fn main() {
    let curve = Curve { a: Fp::new(1), b: Fp::new(1) };
    // 4a^3 + 27b^2 = 4 + 27 = 31, nonzero in F_65521 — curve is non-singular.

    let p1 = find_point(&curve, 0);
    let p2 = find_point(&curve, 5);
    println!("P1 = {:?}", p1);
    println!("P2 = {:?}", p2);
    let (x1, y1) = match p1 { Point::Affine(x, y) => (x, y), _ => panic!() };
    let (x2, y2) = match p2 { Point::Affine(x, y) => (x, y), _ => panic!() };
    assert!(on_curve(&curve, x1, y1));
    assert!(on_curve(&curve, x2, y2));

    // P3 = -(P1 + P2), so P1 + P2 + P3 = O.
    let sum12 = add(&curve, p1, p2);
    let p3 = neg(sum12);
    let (x3, y3) = match p3 { Point::Affine(x, y) => (x, y), _ => panic!() };
    assert!(on_curve(&curve, x3, y3));
    println!("P3 = {:?}", p3);
    println!("P1 + P2 + P3 = {:?}", add(&curve, add(&curve, p1, p2), p3));

    let s2 = semaev_s2(&curve);
    let s3 = semaev_s3(&curve);
    println!("|S_2| = {} terms, |S_3| = {} terms", s2.nterms(), s3.nterms());

    // Sanity: S_2(x1, x1, *) = 0
    let v = s2.eval(&[x1, x1, Fp::zero()]);
    assert!(v.is_zero(), "S_2(x, x, *) should vanish, got {}", v);
    let v = s2.eval(&[x1, x2, Fp::zero()]);
    assert!(!v.is_zero(), "S_2(x1, x2) should be nonzero for x1 != x2");
    println!("S_2 checks OK");

    // Main check: S_3 on (x1, x2, x3) where P1+P2+P3 = O must be zero.
    let v123 = s3.eval(&[x1, x2, x3]);
    let v132 = s3.eval(&[x1, x3, x2]);
    let v231 = s3.eval(&[x2, x3, x1]);
    let v321 = s3.eval(&[x3, x2, x1]);
    println!("S_3(x1,x2,x3) = {}   (expect 0)", v123);
    println!("S_3(x1,x3,x2) = {}   (expect 0, symmetry)", v132);
    println!("S_3(x2,x3,x1) = {}   (expect 0, symmetry)", v231);
    println!("S_3(x3,x2,x1) = {}   (expect 0, symmetry)", v321);
    assert!(v123.is_zero(), "S_3(x1,x2,x3) should vanish");
    assert!(v132.is_zero(), "symmetry: S_3(x1,x3,x2) should vanish");
    assert!(v231.is_zero(), "symmetry: S_3(x2,x3,x1) should vanish");
    assert!(v321.is_zero(), "symmetry: S_3(x3,x2,x1) should vanish");

    // Negative case: a triple that does NOT sum to zero should give nonzero S_3.
    // Use (x1, x2, x1) — P3 = P1, so P1+P2+P1 = 2P1+P2 != O generically.
    let v_bad = s3.eval(&[x1, x2, x1]);
    println!("S_3(x1,x2,x1) = {}   (expect nonzero)", v_bad);
    assert!(!v_bad.is_zero(), "S_3 should be nonzero on non-summing triple");

    println!("\n--- S_4 via resultant ---");
    let s4 = semaev_s4(&curve);
    println!("|S_4| = {} terms, total degree = {}", s4.nterms(), s4.total_degree());

    // Build a true 4-summing tuple (Q1, Q2, Q3, Q4) with Q1+Q2+Q3+Q4 = O.
    // Reusing (p1, p2, p3) would force Q4 = O; pick a fresh third point.
    let q3 = find_point(&curve, 100);
    let q4 = neg(add(&curve, add(&curve, p1, p2), q3));
    let (xq3, _) = match q3 { Point::Affine(x, y) => (x, y), _ => panic!() };
    let (xq4, _) = match q4 { Point::Affine(x, y) => (x, y), _ => panic!("Q4 = O — pick another point") };
    println!("Q1=P1={:?}", p1);
    println!("Q2=P2={:?}", p2);
    println!("Q3   ={:?}", q3);
    println!("Q4   ={:?}", q4);
    println!("Q1+Q2+Q3+Q4 = {:?}",
        add(&curve, add(&curve, add(&curve, p1, p2), q3), q4));

    // S_4 lives in 5-var ambient [X1, X2, X3, X4, Y]. Eval at [x1, x2, xq3, xq4, 0].
    let v4 = s4.eval(&[x1, x2, xq3, xq4, Fp::zero()]);
    println!("S_4(x1,x2,xq3,xq4) = {}   (expect 0)", v4);
    assert!(v4.is_zero(), "S_4 must vanish on 4-summing tuples");

    // Negative case: replace xq4 with x1.
    let v4_bad = s4.eval(&[x1, x2, xq3, x1, Fp::zero()]);
    println!("S_4(x1,x2,xq3,x1)  = {}   (expect nonzero)", v4_bad);
    assert!(!v4_bad.is_zero(), "S_4 should be nonzero off the variety");

    println!("\nAll Semaev S_2 / S_3 / S_4 checks passed.");
}
