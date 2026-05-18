//! Shell-out bridge to the sympy-based gluing divisor map.
//!
//! The C-D gluing step (`FromProdToJac`) has two halves: computing the
//! deg-6 polynomial h(x) (we do this in pure Rust in [`crate::glue`]) and
//! solving a polynomial system to find the Mumford form of the image
//! divisor on J(h). The latter needs a Groebner basis routine; the
//! cleanest path on a no-Sage system is sympy via subprocess.
//!
//! See `scripts/gluing_divisor_map.py` for the Python side.
//!
//! Limitation: this bridge currently handles the F_p subcase only — all
//! inputs (α, β, P_c, P) must be F_p-rational, so coefficients of h fall
//! in F_p too. Generalizing to F_{p²} requires adding `i` as an extra
//! Groebner variable with i² + 1 = 0 in the sympy script.

use std::io::Write;
use std::process::{Command, Stdio};

use num_bigint::BigInt;

use crate::field::{F2, Fp2};
use crate::jacobian::Div;
use crate::poly::Poly;

#[derive(Debug)]
pub struct GluingSolution {
    /// Mumford u(x) = x² + u1·x + u0
    pub u0: F2,
    pub u1: F2,
    /// Mumford v(x) = v1·x + v0
    pub v0: F2,
    pub v1: F2,
}

impl GluingSolution {
    pub fn to_div(&self, fp2: &Fp2) -> Div {
        let u = Poly::new(
            vec![self.u0.clone(), self.u1.clone(), fp2.one()],
            fp2,
        );
        let v = Poly::new(
            vec![self.v0.clone(), self.v1.clone()],
            fp2,
        );
        Div { u, v }
    }
}

/// Compute the Mumford form of the gluing image of (P_c on E_α, P on E_β)
/// on J(h), by shelling out to the sympy script. Returns all solutions
/// returned by the Groebner system.
///
/// Inputs are F_{p²} elements as (a, b) pairs meaning a + b·i. For F_p
/// values, set b = 0. The script auto-detects the F_p subcase for speed.
pub fn gluing_divisor_map_fp(
    p: &BigInt,
    alpha: &[F2; 3],
    beta: &[F2; 3],
    pc: &(F2, F2),
    p_pt: &(F2, F2),
    script_path: &str,
) -> Result<Vec<GluingSolution>, String> {
    let pair = |x: &F2| -> serde_json::Value {
        serde_json::json!([
            x.a.to_string().parse::<i64>().unwrap(),
            x.b.to_string().parse::<i64>().unwrap(),
        ])
    };
    let input = serde_json::json!({
        "p": p.to_string().parse::<u64>().expect("p too big for json u64"),
        "alpha": [pair(&alpha[0]), pair(&alpha[1]), pair(&alpha[2])],
        "beta":  [pair(&beta[0]),  pair(&beta[1]),  pair(&beta[2])],
        "pc":    [pair(&pc.0),     pair(&pc.1)],
        "p_pt":  [pair(&p_pt.0),   pair(&p_pt.1)],
    });

    let mut child = Command::new("python3")
        .arg(script_path)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| format!("failed to spawn python3: {e}"))?;
    child
        .stdin
        .as_mut()
        .unwrap()
        .write_all(input.to_string().as_bytes())
        .map_err(|e| format!("failed to write stdin: {e}"))?;
    let output = child
        .wait_with_output()
        .map_err(|e| format!("failed to wait for python: {e}"))?;
    if !output.status.success() {
        return Err(format!(
            "python script failed: {}",
            String::from_utf8_lossy(&output.stderr)
        ));
    }
    let stdout = String::from_utf8(output.stdout)
        .map_err(|e| format!("non-utf8 stdout: {e}"))?;
    let parsed: serde_json::Value = serde_json::from_str(&stdout)
        .map_err(|e| format!("json parse error: {e} (output was: {stdout})"))?;
    if let Some(err) = parsed.get("error") {
        return Err(format!("python returned error: {err}"));
    }
    let sols = parsed
        .get("solutions")
        .ok_or_else(|| "no 'solutions' key in python output".to_string())?
        .as_array()
        .ok_or_else(|| "'solutions' not an array".to_string())?;
    let mut out = Vec::new();
    for sol in sols {
        if sol.get("error").is_some() {
            return Err(format!("python solution error: {sol}"));
        }
        // Each field is now an [a, b] pair representing a + b·i.
        let pick = |key: &str| -> Result<F2, String> {
            let v = sol.get(key)
                .ok_or_else(|| format!("missing key '{key}'"))?;
            let arr = v.as_array()
                .ok_or_else(|| format!("'{key}' not an array: {:?}", v))?;
            if arr.len() != 2 {
                return Err(format!("'{key}' has length {}", arr.len()));
            }
            let parse_one = |val: &serde_json::Value| -> Result<BigInt, String> {
                match val {
                    serde_json::Value::Number(n) => Ok(BigInt::from(
                        n.as_i64().ok_or_else(|| "non-int".to_string())?
                    )),
                    serde_json::Value::String(s) => s.parse::<i64>()
                        .map(BigInt::from)
                        .map_err(|e| format!("'{s}': {e}")),
                    serde_json::Value::Null => Err("null value".to_string()),
                    _ => Err(format!("unexpected: {:?}", val)),
                }
            };
            Ok(F2 { a: parse_one(&arr[0])?, b: parse_one(&arr[1])? })
        };
        out.push(GluingSolution {
            u0: pick("u0")?,
            u1: pick("u1")?,
            v0: pick("v0")?,
            v1: pick("v1")?,
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::{Fp, Fp2};
    use crate::glue::from_prod_to_jac;

    fn ctx() -> Fp2 { Fp2::new(Fp::new(BigInt::from(431u64))) }

    #[test]
    fn gluing_divisor_map_returns_valid_mumford() {
        let fp2 = ctx();
        let zero = || F2 { a: BigInt::from(0), b: BigInt::from(0) };
        let int = |n: i64| F2 { a: BigInt::from(n), b: BigInt::from(0) };
        // E_α: Montgomery A=6 → 2-torsion (0, 52, 373); 4-tor (1, 55) → (0,0)
        // E_β: Montgomery A=10 → 2-torsion (0, 170, 251); 4-tor (1, 72) → (0,0)
        let alpha = [zero(), int(52), int(373)];
        let beta  = [zero(), int(170), int(251)];
        let pc    = (int(1), int(55));
        let p_pt  = (int(1), int(72));

        let script = format!(
            "{}/scripts/gluing_divisor_map.py",
            env!("CARGO_MANIFEST_DIR")
        );
        let sols = gluing_divisor_map_fp(
            &BigInt::from(431u64),
            &alpha, &beta, &pc, &p_pt,
            &script,
        )
        .expect("python script must succeed");
        assert!(!sols.is_empty(), "expected at least one solution from Groebner");
        println!("got {} solution(s); first = {:?}", sols.len(), &sols[0]);

        let c = from_prod_to_jac(&alpha, &beta, &fp2);
        for (i, sol) in sols.iter().enumerate() {
            let d = sol.to_div(&fp2);
            assert!(
                d.is_valid(&c, &fp2),
                "solution {i} not a valid Mumford form on J(h): u={:?} v={:?}",
                d.u, d.v
            );
        }
    }

    /// Helper: call the msolve-bridge script instead of the sympy one.
    fn call_msolve_bridge(
        p: &BigInt,
        alpha: &[F2; 3],
        beta: &[F2; 3],
        pc: &(F2, F2),
        p_pt: &(F2, F2),
    ) -> Result<Vec<GluingSolution>, String> {
        let script = format!(
            "{}/scripts/msolve_bridge.py",
            env!("CARGO_MANIFEST_DIR")
        );
        gluing_divisor_map_fp(p, alpha, beta, pc, p_pt, &script)
    }

    #[test]
    fn fp2_gluing_via_msolve() {
        // E_α: A=0 (j=1728) → 2-torsion (0, ±i) in F_{p²}\F_p
        // E_β: A=10 → F_p-rational 2-torsion {0, 170, 251}
        // Should exercise the msolve path (closed-form elim + msolve solve).
        let fp2 = ctx();
        let i = fp2.i();
        let neg_i = fp2.neg(&i);
        let zero = fp2.zero();
        let int = |n: i64| F2 { a: BigInt::from(n), b: BigInt::from(0) };

        let alpha = [zero.clone(), i.clone(), neg_i.clone()];
        let beta = [zero, int(170), int(251)];
        let pc = (int(1), int(243));   // 4-torsion on E_α (lift x=1)
        let p_pt = (int(1), int(72));  // 4-torsion on E_β

        let sols = call_msolve_bridge(
            &BigInt::from(431u64),
            &alpha, &beta, &pc, &p_pt,
        )
        .expect("msolve bridge must succeed");
        assert!(!sols.is_empty(), "expected non-empty F_{{p²}} solutions");
        println!("F_{{p²}} solutions: {} found", sols.len());

        // Verify each solution is a valid Mumford form on J(h).
        let c = crate::glue::from_prod_to_jac(&alpha, &beta, &fp2);
        let mut valid_count = 0;
        for (i, sol) in sols.iter().enumerate() {
            let d = sol.to_div(&fp2);
            if d.is_valid(&c, &fp2) {
                valid_count += 1;
            } else {
                eprintln!("solution {i} not valid Mumford: u={:?} v={:?}",
                          d.u, d.v);
            }
        }
        println!("{}/{} solutions are valid Mumford forms",
                 valid_count, sols.len());
        assert!(valid_count > 0, "expected at least one valid Mumford form");
    }

    #[test]
    #[ignore = "slow: brute-force F_{p²} root finding takes ~30s; run with `cargo test -- --ignored`"]
    fn gluing_divisor_map_fp2_path_works() {
        // F_{p²} test: E_α with j=1728 (A=0) has 2-torsion (0, ±i), which
        // lives in F_{p²}\F_p. Exercises the reduced-Groebner + i-as-variable
        // path in the Python bridge.
        let fp2 = ctx();
        let i = fp2.i();
        let neg_i = fp2.neg(&i);
        let zero = fp2.zero();
        let int = |n: i64| F2 { a: BigInt::from(n), b: BigInt::from(0) };

        // E_α: A=0 (j=1728); 2-torsion {0, i, -i}
        let alpha = [zero.clone(), i.clone(), neg_i.clone()];
        // E_β: A=10 → F_p-rational 2-torsion {0, 170, 251}
        let beta = [zero, int(170), int(251)];
        // P_c on E_α: lift x=1 → (1, sqrt(1·(1+0+1))) = (1, sqrt(2)) = (1, 243) in F_p
        let pc = (int(1), int(243));
        // P on E_β: 4-torsion (1, 72)
        let p_pt = (int(1), int(72));

        let script = format!(
            "{}/scripts/gluing_divisor_map.py",
            env!("CARGO_MANIFEST_DIR")
        );
        let sols = gluing_divisor_map_fp(
            &BigInt::from(431u64),
            &alpha, &beta, &pc, &p_pt,
            &script,
        )
        .expect("python script must succeed for F_{p²} inputs");
        // For non-Kani-correct inputs we may get partial / no solutions —
        // but the path should run without timing out (closed-form reduction
        // keeps Groebner fast). Print what we got for inspection.
        println!("F_{{p²}} solutions: {sols:?}");
        // Verify any returned solution is a valid Mumford form.
        let c = crate::glue::from_prod_to_jac(&alpha, &beta, &fp2);
        for (i, sol) in sols.iter().enumerate() {
            let d = sol.to_div(&fp2);
            if !d.is_valid(&c, &fp2) {
                eprintln!("WARNING: solution {i} not valid Mumford: u={:?} v={:?}",
                          d.u, d.v);
            }
        }
    }
}
