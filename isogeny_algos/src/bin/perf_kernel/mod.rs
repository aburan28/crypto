//! Frozen input boundaries for kernel constructors. No constructor algorithm changes.
use super::*;
use isogeny_algos::kernel::{sqrt_velu, velu, xonly};

pub const METHODS: [&str; 7] = [
    "velu_from_point",
    "velu_cyclic",
    "velu_xonly",
    "velu_xonly_fast",
    "sqrt_velu",
    "sqrt_velu_fast",
    "kohel_given_kernel",
];
pub const OUTPUT_CHECK: &str =
    "degree and codomain plus full point transports at 2*ell+1 distinct non-kernel abscissae; Kohel and Velu references";

fn construct<F: Field>(
    f: &F,
    method: &str,
    curve: &Curve<F::E>,
    point: &Pt<F::E>,
    kernel: &Vec<F::E>,
    ell: u64,
) -> Result<Box<dyn Isogeny<F>>, String> {
    let x = match *point {
        Pt::Aff(x, _) => x,
        Pt::Inf => return Err("infinite kernel generator".into()),
    };
    Ok(match method {
        "velu_from_point" => Box::new(velu::velu_from_point(f, curve, point, ell)),
        "velu_cyclic" => Box::new(velu::velu_cyclic(f, curve, point, ell)),
        "velu_xonly" => {
            Box::new(xonly::velu_xonly(f, curve, x, ell).ok_or("x-only constructor failed")?)
        }
        "velu_xonly_fast" => Box::new(
            xonly::velu_xonly_fast(f, curve, x, ell).ok_or("fast x-only constructor failed")?,
        ),
        "sqrt_velu" => Box::new(sqrt_velu::sqrt_velu(f, curve, point, ell)),
        "sqrt_velu_fast" => Box::new(
            sqrt_velu::sqrt_velu_fast(f, curve, x, ell)
                .ok_or("fast square-root constructor failed")?,
        ),
        "kohel_given_kernel" => Box::new(kohel(f, curve, kernel, ell)),
        _ => return Err("unknown kernel constructor".into()),
    })
}
fn point_record<const N: usize>(f: &FpM<N>, p: &Pt<[u64; N]>) -> Json {
    match p {
        Pt::Inf => Json::Null,
        Pt::Aff(x, y) => Json::Arr(vec![
            Json::str(f.to_big(x).to_dec()),
            Json::str(f.to_big(y).to_dec()),
        ]),
    }
}
fn output<F: Field<E = [u64; N]>, const N: usize>(
    canonical: &FpM<N>,
    f: &F,
    iso: &dyn Isogeny<F>,
    probes: &[Pt<[u64; N]>],
) -> Json {
    Json::obj(vec![
        ("degree", Json::str(iso.degree())),
        ("a", Json::str(canonical.to_big(&iso.codomain().a).to_dec())),
        ("b", Json::str(canonical.to_big(&iso.codomain().b).to_dec())),
        (
            "point_transports",
            Json::Arr(
                probes
                    .iter()
                    .map(|p| point_record(canonical, &iso.eval(f, p)))
                    .collect(),
            ),
        ),
    ])
}
pub fn run<const N: usize>(
    out: &mut std::fs::File,
    bits: usize,
    ells: &[u64],
    samples: usize,
    timing: bool,
) -> Result<(), String> {
    let mut primes = vec![3];
    primes.extend_from_slice(ells);
    primes.sort();
    primes.dedup();
    let mut rng = Rng::new(0x4b45_524e_454c + bits as u64);
    let w = supersingular_workload::<N>(&primes, bits, &mut rng);
    for &ell in ells {
        let p = w
            .pts
            .iter()
            .find(|(l, _)| *l == ell)
            .ok_or("missing generator")?
            .1;
        let (g, _) = kernel_poly_from_point_f(&w.f, &w.e, &p, ell);
        let reference = kohel(&w.f, &w.e, &g, ell);
        let mut probes = vec![];
        let mut seen = std::collections::BTreeSet::new();
        while probes.len() < (2 * ell + 1) as usize {
            let q = random_point_f(&w.f, &w.e, &mut rng);
            if let Pt::Aff(x, _) = q {
                if reference.eval(&w.f, &q) != Pt::Inf && seen.insert(x) {
                    probes.push(q);
                }
            }
        }
        let expected = output(&w.f, &w.f, &reference, &probes).dump();
        let point_reference = velu::velu_from_point(&w.f, &w.e, &p, ell);
        if output(&w.f, &w.f, &point_reference, &probes).dump() != expected {
            return Err("independent Kohel and Velu references disagree".into());
        }
        for method in METHODS {
            eprintln!("kernel {bits}-bit ell={ell} {method}");
            let boundary = match method {
                "kohel_given_kernel" => "supplied kernel polynomial -> owned rational map; deriving kernel excluded",
                "velu_xonly" | "velu_xonly_fast" | "sqrt_velu_fast" => "supplied generator x -> owned evaluation representation; multiples included",
                _ => "supplied affine generator -> owned evaluation representation; multiples included",
            };
            let input = Json::obj(vec![
                ("p", Json::str(w.f.q().to_dec())),
                ("a", Json::str(w.f.to_big(&w.e.a).to_dec())),
                ("b", Json::str(w.f.to_big(&w.e.b).to_dec())),
                ("ell", Json::str(ell)),
                ("generator", point_record(&w.f, &p)),
                ("reference_kernel", elements(&w.f, &g)),
                ("reference_map_sha256", Json::str(sha256_hex(map_record(&w.f, &reference).dump().as_bytes()))),
                ("probes", Json::Arr(probes.iter().map(|q| point_record(&w.f, q)).collect())),
                ("boundary", Json::str(boundary)),
                ("excluded", Json::str("field/fixture generation, reference, fingerprints and verification; result allocation included in constructor")),
            ]);
            let counted = Counted::new(w.f.clone());
            let counted_iso = construct(&counted, method, &w.e, &p, &g, ell)?;
            // Save counts before transport verification touches the counted field.
            let counts = counted.snapshot();
            if output(&w.f, &counted, counted_iso.as_ref(), &probes).dump() != expected {
                return Err(format!("{method} counted output mismatch"));
            }
            let native = construct(&w.f, method, &w.e, &p, &g, ell)?;
            if output(&w.f, &w.f, native.as_ref(), &probes).dump() != expected {
                return Err(format!("{method} native output mismatch"));
            }
            let mut raw = vec![];
            if timing {
                black_box(construct(&w.f, method, &w.e, &p, &g, ell)?);
                for _ in 0..samples {
                    let started = Instant::now();
                    let iso = black_box(construct(
                        black_box(&w.f),
                        method,
                        black_box(&w.e),
                        black_box(&p),
                        black_box(&g),
                        ell,
                    )?);
                    let ns = started.elapsed().as_nanos();
                    if output(&w.f, &w.f, iso.as_ref(), &probes).dump() != expected {
                        return Err(format!("{method} timed output mismatch"));
                    }
                    raw.push(Json::str(ns));
                }
            }
            let row = Json::obj(vec![
                ("schema", Json::str("isogeny-perf.case/v1")),
                ("bits", Json::str(bits)),
                ("ell", Json::str(ell)),
                ("algorithm", Json::str(method)),
                (
                    "input_sha256",
                    Json::str(sha256_hex(input.dump().as_bytes())),
                ),
                ("input", input),
                ("output_sha256", Json::str(sha256_hex(expected.as_bytes()))),
                ("output_check", Json::str(OUTPUT_CHECK)),
                ("count_model", Json::str(COUNT_MODEL)),
                ("operations", counts),
                ("diagnostic_ns", Json::Arr(raw)),
                ("qualified_timing", Json::Bool(false)),
                ("verified", Json::Bool(true)),
            ]);
            writeln!(out, "{}", row.dump()).map_err(|e| e.to_string())?;
            out.flush().map_err(|e| e.to_string())?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn all_seven_constructors_match_an_independent_map_on_public_points() {
        let f = isogeny_algos::field::Zp::new(101);
        let mut rng = Rng::new(451);
        let (e, p) = isogeny_algos::testdata::curve_with_point(&f, 7, &mut rng);
        let (kernel, _) = kernel_poly_from_point(&f, &e, &p, 7);
        let reference = kohel(&f, &e, &kernel, 7);
        let points: Vec<_> = (0..20).map(|_| e.random_point(&f, &mut rng)).collect();
        for method in METHODS {
            let iso = construct(&f, method, &e, &p, &kernel, 7).unwrap();
            assert_eq!(iso.degree(), 7, "{method}");
            assert_eq!(iso.codomain(), &reference.cod, "{method}");
            for q in &points {
                assert_eq!(iso.eval(&f, q), reference.eval(&f, q), "{method}");
            }
        }
    }
}
