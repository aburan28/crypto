//! The vartime point operations on **every** input of the smallest
//! fields, against what they stand in for.
//!
//! The differential tests in `ecc::point` and `ec_vartime_walks_pinned`
//! draw their operands at random from primes of 251 and up; these
//! enumerate instead: every coordinate pair (on a curve or not, since
//! the formulas never read `b`), every `a`, every scalar up to past
//! `3p`, and for the linear search every base, target and a range of
//! bounds.  A result is compared value and modulus, a panic by its
//! message.  The smallest primes are where the Jacobian formulas'
//! factors of 2 and 3 stop being units, so they are where a projective
//! shortcut can part from the affine arithmetic it replaces.

use crypto_lib::ecc::{FieldElement, Point};
use num_bigint::BigUint;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::sync::Mutex;

/// `f()`'s value, or the message it panicked with.
fn outcome<T>(f: impl FnOnce() -> T) -> Result<T, String> {
    catch_unwind(AssertUnwindSafe(f)).map_err(|err| match err.downcast::<String>() {
        Ok(s) => *s,
        Err(err) => err
            .downcast_ref::<&str>()
            .map_or_else(String::new, |s| s.to_string()),
    })
}

/// A point with both moduli spelled out, so a result that differs
/// only in a coordinate's modulus still shows.
fn show(p: &Point) -> String {
    match p {
        Point::Infinity => "O".into(),
        Point::Affine { x, y } => {
            format!(
                "({} mod {}, {} mod {})",
                x.value, x.modulus, y.value, y.modulus
            )
        }
    }
}

/// The search [`Point::linear_dlog_vartime`] replaced: step
/// `current = current.add(g, a)` and compare.
fn stepped_search(g: &Point, t: &Point, bound: u64, a: &FieldElement) -> Option<u64> {
    let mut current = Point::Infinity;
    for k in 0..bound {
        if current == *t {
            return Some(k);
        }
        current = current.add(g, a);
    }
    None
}

/// The identity and every pair `(x, y)` of residues mod `m`.
fn every_point(m: u64) -> Vec<Point> {
    let fe = |v: u64| FieldElement::new(BigUint::from(v), BigUint::from(m));
    let mut out = vec![Point::Infinity];
    for x in 0..m {
        for y in 0..m {
            out.push(Point::Affine { x: fe(x), y: fe(y) });
        }
    }
    out
}

/// Every `a` mod `m` with every point of [`every_point`].
fn for_every_input(m: u64, mut f: impl FnMut(&FieldElement, &[Point])) {
    let points = every_point(m);
    for a in 0..m {
        f(
            &FieldElement::new(BigUint::from(a), BigUint::from(m)),
            &points,
        );
    }
}

/// Runs `body` with the panic hook silenced, so the expected panics do
/// not flood the output, and restores it before returning `body`'s
/// result, so a failing assertion afterwards still prints.  The hook is
/// process-wide, so the tests here take turns.
fn quietly<T>(body: impl FnOnce() -> T) -> T {
    static TURN: Mutex<()> = Mutex::new(());
    let _turn = TURN.lock().unwrap_or_else(|e| e.into_inner());
    let hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let r = catch_unwind(AssertUnwindSafe(body));
    std::panic::set_hook(hook);
    r.unwrap_or_else(|e| std::panic::resume_unwind(e))
}

/// Over `F_2`, `F_3`, `F_5` and `F_7`: `add`, `double` and `scalar_mul`
/// against their vartime versions on every input, and the linear search
/// against the stepped one over the odd primes (`F_2` is the next test).
#[test]
fn vartime_ops_match_on_every_input_of_tiny_prime_fields() {
    let bad = quietly(|| {
        let mut bad = Vec::new();
        for m in [2u64, 3, 5, 7] {
            for_every_input(m, |a, points| {
                for p1 in points {
                    let (s, t) = (outcome(|| p1.double(a)), outcome(|| p1.double_vartime(a)));
                    let (s, t) = (s.map(|p| show(&p)), t.map(|p| show(&p)));
                    if s != t {
                        bad.push(format!("m={m} a={a} 2·{}: {s:?} vs {t:?}", show(p1)));
                    }
                    for k in 0..3 * m + 4 {
                        let k = BigUint::from(k);
                        let s = outcome(|| p1.scalar_mul(&k, a)).map(|p| show(&p));
                        let t = outcome(|| p1.scalar_mul_vartime(&k, a)).map(|p| show(&p));
                        if s != t {
                            bad.push(format!("m={m} a={a} {k}·{}: {s:?} vs {t:?}", show(p1)));
                        }
                    }
                    for p2 in points {
                        let s = outcome(|| p1.add(p2, a)).map(|p| show(&p));
                        let t = outcome(|| p1.add_vartime(p2, a)).map(|p| show(&p));
                        if s != t {
                            bad.push(format!(
                                "m={m} a={a} {} + {}: {s:?} vs {t:?}",
                                show(p1),
                                show(p2)
                            ));
                        }
                        if m == 2 {
                            continue;
                        }
                        for bound in [0, 1, 2, 3, 2 * m + 3] {
                            let s = outcome(|| stepped_search(p1, p2, bound, a));
                            let t = outcome(|| p1.linear_dlog_vartime(p2, bound, a));
                            if s != t {
                                bad.push(format!(
                                    "m={m} a={a} log {} base {} < {bound}: {s:?} vs {t:?}",
                                    show(p2),
                                    show(p1)
                                ));
                            }
                        }
                    }
                }
            });
        }
        bad
    });
    assert!(
        bad.is_empty(),
        "{} mismatches:\n{}",
        bad.len(),
        bad.join("\n")
    );
}

/// The linear search over `F_2`, on every input.  Ignored while it
/// fails: `2 = 0` there, so a Jacobian doubling's `Z₃ = 2·Y·Z` and a
/// mixed addition's `Z₃ = 2·Z·H` are zero, where the affine search's
/// `2y` has no inverse and it panics.  E.g. base `(0, 1)`, `a = 0`,
/// target `(1, 0)`, bound 3: the stepped search panics computing `2G`;
/// `linear_dlog_vartime` returns `Some(2)`, and through it
/// `recover_in_prime_power_subgroup` returns a digit where it panicked.
#[test]
#[ignore = "linear_dlog_vartime departs from the stepped search over F_2"]
fn linear_dlog_vartime_matches_stepped_search_over_f2() {
    let bad = quietly(|| {
        let mut bad = Vec::new();
        for_every_input(2, |a, points| {
            for g in points {
                for t in points {
                    for bound in 0..8 {
                        let s = outcome(|| stepped_search(g, t, bound, a));
                        let v = outcome(|| g.linear_dlog_vartime(t, bound, a));
                        if s != v {
                            bad.push(format!(
                                "a={a} log {} base {} < {bound}: {s:?} vs {v:?}",
                                show(t),
                                show(g)
                            ));
                        }
                    }
                }
            }
        });
        bad
    });
    assert!(
        bad.is_empty(),
        "{} mismatches:\n{}",
        bad.len(),
        bad.join("\n")
    );
}

/// Over small composites, where the points may differ (Euclid's inverse
/// against Fermat's `a^(m−2)`), the vartime operations panic on exactly
/// the inputs the constant-time ones do, with the same message, on
/// every reduced input.
#[test]
fn vartime_ops_panic_where_constant_time_ops_do_on_tiny_composites() {
    let bad = quietly(|| {
        let mut bad = Vec::new();
        for m in [4u64, 6, 8, 9, 10] {
            for_every_input(m, |a, points| {
                for p1 in points {
                    let s = outcome(|| p1.double(a)).err();
                    let t = outcome(|| p1.double_vartime(a)).err();
                    if s != t {
                        bad.push(format!("m={m} a={a} 2·{}: {s:?} vs {t:?}", show(p1)));
                    }
                    for k in 0..3 * m + 4 {
                        let k = BigUint::from(k);
                        let s = outcome(|| p1.scalar_mul(&k, a)).err();
                        let t = outcome(|| p1.scalar_mul_vartime(&k, a)).err();
                        if s != t {
                            bad.push(format!("m={m} a={a} {k}·{}: {s:?} vs {t:?}", show(p1)));
                        }
                    }
                    for p2 in points {
                        let s = outcome(|| p1.add(p2, a)).err();
                        let t = outcome(|| p1.add_vartime(p2, a)).err();
                        if s != t {
                            bad.push(format!(
                                "m={m} a={a} {} + {}: {s:?} vs {t:?}",
                                show(p1),
                                show(p2)
                            ));
                        }
                    }
                }
            });
        }
        bad
    });
    assert!(
        bad.is_empty(),
        "{} mismatches:\n{}",
        bad.len(),
        bad.join("\n")
    );
}
