//! Executable prerequisites, kept separate from measured ecbench pipelines.
use crypto_lib::binary_ecc::{
    curve::scalar_mul, BinaryCurve, BinaryPoint, F2mElement, MumfordDivisor, OrdinaryBinaryCover,
};
use crypto_lib::cryptanalysis::{
    ecbench::{
        conductor_gap,
        workload::{registry_koblitz_modulus, CurveSpec, Instance},
    },
    ic_boundary::{koblitz_instance, BinaryGroup, CountedGroup, GroupOps},
    ic_framework::{
        plugins::CyclotomicBase,
        stages::{FactorBaseBuilder, InstanceCtx, Params},
    },
    isogeny_degree_search::{full_field_basis, present, rank_profile},
    koblitz_groebner::FieldStructure,
    koblitz_index_calculus::find_irreducible,
    koblitz_wide::WideKoblitz,
    pq_groebner_f2::F2BoolPoly,
    quasi_subfield::{census, cyclotomic_cosets, stable_subspace_dimensions},
    semaev_decomp::Gf2,
};
use num_bigint::BigUint;
use serde_json::{json, Value};
use std::collections::BTreeMap;

fn err(e: impl std::fmt::Display) -> String {
    e.to_string()
}

fn quasi() -> Result<Value, String> {
    let mut rows = Vec::new();
    for (n, n0) in [(5, 2), (7, 3), (7, 4), (11, 3)] {
        let c = census(n, n0, 3).ok_or("census refused frozen input")?;
        let irr = find_irreducible(n).ok_or("no irreducible polynomial")?;
        let gf = Gf2::new(&irr);
        let curve = koblitz_instance(0, n);
        let mut witnesses = Vec::new();
        for w in &c.witnesses {
            let roots: Vec<u64> = (0..(1u64 << n))
                .filter(|&x| {
                    let mut y = gf.sqr_k(x, n0);
                    for (i, &a) in w.lambda.iter().enumerate() {
                        y ^= gf.mul(a, gf.sqr_k(x, i as u32));
                    }
                    y == 0
                })
                .collect();
            if roots.len() != 1usize << n0 {
                return Err("rank/root-count disagreement".into());
            }
            let points: Vec<_> = match &curve {
                Some(inst) => roots.iter().flat_map(|&x| inst.points_with_x(x)).collect(),
                None => Vec::new(),
            };
            let in_subgroup = curve.as_ref().map(|inst| {
                points
                    .iter()
                    .filter(|&&p| inst.fast.mul_u64(p, inst.r).infinity)
                    .count()
            });
            witnesses.push(
                json!({"lambda":w.lambda,"j":w.j,"roots":roots,"root_count":roots.len(),
                "curve":curve.as_ref().map(|i|i.name.clone()),"curve_points":points.len(),
                "unprojected_subgroup_points":in_subgroup,"independent_root_check":true}),
            );
        }
        rows.push(json!({"n":n,"n0":n0,"modulus_low_terms":irr.low_terms,"j_max":c.j_max,
            "coefficient_vectors_examined":c.examined,"split_count":c.found,"min_j":c.min_j,
            "witnesses":witnesses,"scope":"complete finite coefficient search, not a scalable family or full IC measurement"}));
    }
    Ok(json!(rows))
}

fn nonlinear() -> Result<Value, String> {
    let mut rows = Vec::new();
    for n in [11, 23] {
        let inst = koblitz_instance(1, n).ok_or("missing nonlinear instance")?;
        let g = BinaryGroup(&inst.fast);
        let ctx = InstanceCtx {
            group: &g,
            generator: inst.generator,
            target: inst.generator,
            r: inst.r,
            cofactor: inst.cofactor,
            group_order: inst.group_order,
            name: inst.name.clone(),
            field_degree: Some(n),
        };
        for w in [1, 2, 3] {
            let params = Params(BTreeMap::from([
                ("w".into(), w.to_string()),
                ("columns".into(), "8".into()),
            ]));
            let mut ops = GroupOps::default();
            match (CyclotomicBase { instance: &inst }).build(&ctx, &params, &mut ops) {
                Ok(fb) => {
                    let mut all_on = true;
                    let mut all_order = true;
                    let mut closed = true;
                    let mut coefficients = true;
                    let lambda = inst.koblitz.as_ref().unwrap().lambda.to_u64_digits()[0];
                    for (idx, &p) in fb.points.iter().enumerate() {
                        all_on &= inst.fast.is_on_curve(p);
                        all_order &= inst.fast.mul_u64(p, inst.r).infinity;
                        let image = inst.fast.frobenius_k(p, 1);
                        if let Some(j) = fb.index_of_key(g.key(&image)) {
                            coefficients &= fb.col_of[idx] == fb.col_of[j]
                                && ((fb.coef_of[idx] as u128 * lambda as u128) % inst.r as u128)
                                    as u64
                                    == fb.coef_of[j];
                        } else {
                            closed = false;
                        }
                    }
                    if !all_on || !all_order || !closed || !coefficients {
                        return Err("nonlinear base witness failed".into());
                    }
                    rows.push(json!({"curve":inst.name,"n":n,"gamma_weight":w,"status":"verified",
                        "columns":fb.columns,"signed_points":fb.points.len(),"abscissae":fb.abscissae,
                        "construction_ops":{"adds":ops.adds,"doubles":ops.doubles},"counters":fb.cost.native,
                        "on_curve":all_on,"subgroup_annihilation":all_order,"frobenius_closed":closed,
                        "orbit_coefficient_witnesses":coefficients,"normal_basis_hamming_ideal":false}));
                }
                Err(e) => rows.push(
                    json!({"curve":inst.name,"n":n,"gamma_weight":w,"status":"error","error":e}),
                ),
            }
        }
    }
    Ok(json!(rows))
}

fn degree_profiles() -> Result<Value, String> {
    let mut rows = Vec::new();
    for n in [13, 19, 23, 31] {
        let inst = koblitz_instance(0, n).ok_or("missing degree instance")?;
        let st = FieldStructure::new(n, &inst.irreducible);
        for ell in [3, 4] {
            let basis = &full_field_basis(n)[..ell];
            for k in [1, 3] {
                let p = inst.fast.mul_u64(inst.generator, k);
                let xr = F2mElement::from_biguint(&BigUint::from(p.x), n);
                let sys = present(basis, &xr, &F2mElement::one(n), 2, &st)
                    .ok_or("system construction refused")?;
                let chain = present(basis, &xr, &F2mElement::one(n), 3, &st)
                    .ok_or("chained system construction refused")?;
                let profile = rank_profile(&sys, 5);
                let roots = (0..1u64 << sys.n_vars)
                    .filter(|&bits| sys.equations.iter().all(|p| !eval(p, bits)))
                    .count();
                rows.push(json!({"curve":inst.name,"n":n,"ell":ell,"summands":2,"vars":sys.n_vars,
                    "three_summand_chain_variables":chain.n_vars,"three_summand_chain_equations":chain.equations.len(),
                    "three_summand_chain_equations_sha256":crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex(&serde_json::to_vec(&chain.equations).unwrap()),
                    "structural_target_multiple":k,"first_fall_degree_up_to_5":profile.fall_degree,
                    "exhaustive_boolean_roots":roots,"assignments_examined":1u64<<sys.n_vars,
                    "profiles":profile.profiles.iter().map(|p|json!({"degree":p.degree,"rows":p.rows,"columns":p.cols,"rank":p.rank,"row_syzygies":p.syzygies()})).collect::<Vec<_>>(),
                    "scope":"operational Macaulay ranks; neither solving degree nor lifted point-decomposition success"}));
            }
        }
    }
    Ok(json!(rows))
}

fn eval(poly: &F2BoolPoly, bits: u64) -> bool {
    poly.terms
        .iter()
        .fold(false, |v, m| v ^ (bits & m.mask == m.mask))
}

#[derive(Default)]
struct Qubo {
    vars: usize,
    offset: i64,
    coefficients: BTreeMap<(usize, usize), i64>,
    gates: Vec<(usize, usize, usize)>,
    parity: Vec<(i64, Vec<usize>, Vec<usize>)>,
}

impl Qubo {
    fn term(&mut self, i: usize, j: usize, c: i64) {
        *self.coefficients.entry((i.min(j), i.max(j))).or_default() += c;
    }
    fn var(&mut self) -> usize {
        let i = self.vars;
        self.vars += 1;
        i
    }
    fn and(&mut self, x: usize, y: usize) -> usize {
        let z = self.var();
        self.term(x, y, 1);
        self.term(x, z, -2);
        self.term(y, z, -2);
        self.term(z, z, 3);
        self.gates.push((x, y, z));
        z
    }
    fn square(&mut self, c: i64, terms: &[(usize, i64)]) {
        self.offset += c * c;
        for (k, &(i, a)) in terms.iter().enumerate() {
            self.term(i, i, a * a + 2 * c * a);
            for &(j, b) in &terms[..k] {
                self.term(i, j, 2 * a * b);
            }
        }
    }
    fn encode(eqs: &[F2BoolPoly], vars: usize) -> Self {
        let mut q = Self {
            vars,
            ..Self::default()
        };
        let mut products = BTreeMap::new();
        for i in 0..vars {
            products.insert(1u64 << i, i);
        }
        for eq in eqs {
            let mut outputs = Vec::new();
            let mut constant = 0;
            for monomial in &eq.terms {
                if monomial.mask == 0 {
                    constant = 1;
                    continue;
                }
                let mut mask = 0;
                let mut output = None;
                for i in 0..vars {
                    if monomial.mask >> i & 1 == 0 {
                        continue;
                    }
                    mask |= 1u64 << i;
                    let index = if let Some(&v) = products.get(&mask) {
                        v
                    } else {
                        let v = q.and(output.unwrap(), i);
                        products.insert(mask, v);
                        v
                    };
                    output = Some(index);
                }
                outputs.push(output.unwrap());
            }
            let max_slack = (outputs.len() + constant as usize) / 2;
            let mut slack = Vec::new();
            let mut power = 1;
            while power <= max_slack {
                slack.push(q.var());
                power *= 2;
            }
            let mut terms: Vec<_> = outputs.iter().map(|&v| (v, 1i64)).collect();
            terms.extend(
                slack
                    .iter()
                    .enumerate()
                    .map(|(i, &v)| (v, -2 * (1i64 << i))),
            );
            q.square(constant, &terms);
            q.parity.push((constant, outputs, slack));
        }
        q.coefficients.retain(|_, c| *c != 0);
        q
    }
    fn energy(&self, bits: &[bool]) -> i64 {
        self.offset
            + self
                .coefficients
                .iter()
                .filter(|((i, j), _)| bits[*i] && bits[*j])
                .map(|(_, c)| *c)
                .sum::<i64>()
    }
    fn faithful_completion(&self, original: u64, n: usize) -> Vec<bool> {
        let mut bits = vec![false; self.vars];
        for (i, b) in bits.iter_mut().take(n).enumerate() {
            *b = original >> i & 1 == 1;
        }
        for &(x, y, z) in &self.gates {
            bits[z] = bits[x] && bits[y];
        }
        for (c, outputs, slack) in &self.parity {
            let k = (*c as usize + outputs.iter().filter(|&&v| bits[v]).count()) / 2;
            for (i, &v) in slack.iter().enumerate() {
                bits[v] = k >> i & 1 == 1;
            }
        }
        bits
    }
}

fn qubo_control() -> Result<Value, String> {
    let inst = koblitz_instance(0, 7).ok_or("missing QUBO curve")?;
    let st = FieldStructure::new(7, &inst.irreducible);
    let sys = present(
        &full_field_basis(7)[..2],
        &F2mElement::from_biguint(&BigUint::from(inst.generator.x), 7),
        &F2mElement::one(7),
        3,
        &st,
    )
    .ok_or("QUBO system construction failed")?;
    let q = Qubo::encode(&sys.equations, sys.n_vars);
    let mut zero = 0;
    for bits in 0..1u64 << sys.n_vars {
        let completion = q.faithful_completion(bits, sys.n_vars);
        let expected = sys.equations.iter().filter(|p| eval(p, bits)).count() as i64;
        let energy = q.energy(&completion);
        if energy != expected {
            return Err(format!(
                "QUBO energy mismatch at {bits}: {energy} != {expected}"
            ));
        }
        if energy == 0 {
            zero += 1;
        }
    }
    let gates_verified = and_table();
    if !gates_verified {
        return Err("Rosenberg AND truth table failed".into());
    }
    Ok(
        json!({"curve":inst.name,"original_variables":sys.n_vars,"equations":sys.equations,
        "total_variables":q.vars,"auxiliary_variables":q.vars-sys.n_vars,"and_gates":q.gates.len(),
        "coefficient_count":q.coefficients.len(),"offset":q.offset,
        "coefficients":q.coefficients.iter().map(|(&(i,j),&c)|json!([i,j,c])).collect::<Vec<_>>(),
        "original_assignments_examined":1u64<<sys.n_vars,"zero_energy_faithful_completions":zero,
        "faithful_completion_matches_violated_equation_count":true,"and_truth_table_verified":gates_verified,
        "all_auxiliary_assignments_enumerated":false,"hardware_annealing":false,
        "scope":"finite encoding control only; no embedding, sampling or DLP speed claim"}),
    )
}

fn and_table() -> bool {
    (0..8).all(|v| {
        let x = v & 1;
        let y = v >> 1 & 1;
        let z = v >> 2 & 1;
        let p = x * y - 2 * x * z - 2 * y * z + 3 * z;
        p >= 0 && (p == 0) == (z == x * y)
    })
}

fn covers() -> Result<Value, String> {
    let mut rows = Vec::new();
    for (n, a) in [(7, 0), (83, 0), (83, 1)] {
        let (irr, r, h, generator, slug) = if n == 7 {
            let inst = koblitz_instance(a, n).ok_or("missing cover instance")?;
            (
                inst.irreducible,
                BigUint::from(inst.r),
                BigUint::from(inst.cofactor),
                (
                    BigUint::from(inst.generator.x),
                    BigUint::from(inst.generator.y),
                ),
                inst.name,
            )
        } else if a == 1 {
            // E1 at degree83 has a composite #E/2; use the exact catalogued
            // prime subgroup, never silently assert cofactor two.
            let catalog: Value =
                serde_json::from_str(include_str!("../../docs/curves/registry.json"))
                    .map_err(err)?;
            let row = catalog["curves"]
                .as_array()
                .ok_or("missing catalog curves")?
                .iter()
                .find(|c| {
                    c["family"] == "koblitz"
                        && c["params"]["n"] == 83
                        && c["params"]["a"] == 1
                        && c["representations"]
                            .as_array()
                            .is_some_and(|r| !r.is_empty())
                })
                .ok_or("missing exact E1 degree83 subgroup representation")?;
            let rep = &row["representations"][0]["curve"];
            let parse = |s: &Value, radix| {
                BigUint::parse_bytes(
                    s.as_str()
                        .ok_or("expected integer string")?
                        .trim_start_matches("0x")
                        .as_bytes(),
                    radix,
                )
                .ok_or("invalid integer".to_string())
            };
            (
                registry_koblitz_modulus(a, n).ok_or("missing E1 modulus")?,
                parse(&rep["subgroup_order"], 10)?,
                BigUint::from(rep["cofactor"].as_u64().ok_or("missing cofactor")?),
                (
                    parse(&rep["generator"][0], 16)?,
                    parse(&rep["generator"][1], 16)?,
                ),
                row["slug"].as_str().ok_or("missing slug")?.to_string(),
            )
        } else {
            let irr = registry_koblitz_modulus(a, n).ok_or("unregistered wide cover")?;
            let inst = WideKoblitz::new(a, &irr)?;
            let curve = CurveSpec::Koblitz { a, n };
            let Instance::Wide(built) = curve.build()? else {
                return Err("expected wide instance".into());
            };
            let slug = Instance::Wide(built).facts(&curve).slug;
            (
                irr,
                BigUint::from(inst.r),
                BigUint::from(inst.cofactor),
                (
                    BigUint::from(inst.generator.0),
                    BigUint::from(inst.generator.1),
                ),
                slug,
            )
        };
        let g = BinaryPoint::Affine {
            x: F2mElement::from_biguint(&generator.0, n),
            y: F2mElement::from_biguint(&generator.1, n),
        };
        let curve = BinaryCurve {
            m: n,
            irreducible: irr.clone(),
            a: F2mElement::from_biguint(&BigUint::from(a), n),
            b: F2mElement::one(n),
            generator: g.clone(),
            order: r.clone(),
            cofactor: h,
        };
        let cover = OrdinaryBinaryCover::new(n, irr.clone(), curve.a.clone(), curve.b.clone())
            .map_err(err)?;
        cover.verify_certificate().map_err(err)?;
        let source_annihilation = scalar_mul(&curve, &g, &r) == BinaryPoint::Infinity;
        let image = cover.pullback_point_class(&g).map_err(err)?;
        let image_annihilation = image
            .try_scalar_mul(&r, cover.source_curve())
            .map_err(err)?
            .is_identity();
        let q = scalar_mul(&curve, &g, &BigUint::from(7u32));
        let image_q = cover.pullback_point_class(&q).map_err(err)?;
        let scalar_homomorphism = image
            .try_scalar_mul(&BigUint::from(7u32), cover.source_curve())
            .map_err(err)?
            == image_q;
        let composition = cover.pull_push_composition_certificate(&g).map_err(err)?;
        let infinity = cover
            .pullback_point_class(&BinaryPoint::Infinity)
            .map_err(err)?
            .is_identity();
        let torsion = BinaryPoint::Affine {
            x: F2mElement::zero(n),
            y: F2mElement::one(n),
        };
        let torsion_check = cover
            .pull_push_composition_certificate(&torsion)
            .map_err(err)?
            .identity_holds;
        let off_curve = BinaryPoint::Affine {
            x: F2mElement::zero(n),
            y: F2mElement::zero(n),
        };
        let off_curve_rejected = cover.pullback_point_class(&off_curve).is_err();
        let no_general_norm = cover
            .pushforward_arbitrary_jacobian_class(&MumfordDivisor::identity(3))
            .is_err();
        if !source_annihilation
            || image.is_identity()
            || !image_annihilation
            || !scalar_homomorphism
            || !composition.identity_holds
            || !infinity
            || !torsion_check
            || !off_curve_rejected
            || !no_general_norm
        {
            return Err(format!("cover correctness control failed on {slug}"));
        }
        rows.push(json!({"curve":slug,"n":n,"a":a,"modulus_low_terms":irr.low_terms,
            "subgroup_order":r.to_string(),"generator":[generator.0.to_string(),generator.1.to_string()],
            "certificate":format!("{:?}",cover.certificate()),"source_generator_nonidentity":true,"source_order_annihilation":source_annihilation,
            "image_nonidentity":!image.is_identity(),"image_order_annihilation":image_annihilation,
            "sample_scalar":7,"sample_scalar_homomorphism":scalar_homomorphism,"typed_push_pull_identity":composition.identity_holds,
            "infinity_check":infinity,"ramified_two_torsion_check":torsion_check,"off_curve_rejected":off_curve_rejected,
            "general_norm_unsupported":no_general_norm,"subgroup_status":format!("{:?}",cover.subgroup_transfer_status(&r).map_err(err)?),
            "independent_exact_prime_order_certificate":false,"destination_dlp_run":false,
            "scope":"actual genus-3 divisor pullbacks, not a higher-genus index-calculus speed result"}));
    }
    Ok(json!(rows))
}

fn isogeny_controls() -> Result<Value, String> {
    let p = conductor_gap::params(263, 241)?;
    let crater = conductor_gap::crater(p)?;
    let floor = conductor_gap::floor(p, 59297)?;
    let discriminant = 241i128 * 241 - 4 * p.q as i128;
    let valid = discriminant == -19 * 263i128 * 263;
    let annihilates = |i: &crypto_lib::cryptanalysis::ic_boundary::PrimeInstance| {
        let g = crypto_lib::cryptanalysis::ic_boundary::PrimePoint::affine(
            i.generator.0,
            i.generator.1,
        );
        i.curve.mul(&mut GroupOps::default(), g, i.r).infinity
    };
    if !valid || !annihilates(&crater) || !annihilates(&floor) {
        return Err("CM controls failed".into());
    }
    Ok(
        json!({"ell":263,"trace":241,"field_order":p.q,"subgroup_order":p.n,
        "crater":crater.name,"floor":floor.name,"discriminant":discriminant,
        "discriminant_identity":valid,"source_and_destination_generators_annihilated":true,
        "bad_even_trace_rejected":conductor_gap::params(263,242).is_err(),
        "bad_prime_pair_rejected":conductor_gap::params(263,3).is_err(),
        "explicit_isogeny_map":false,"cross_model_target_pairing":false,"hecke_rank_theorem":false}),
    )
}

pub fn run() -> Result<Value, String> {
    let mut module = Vec::new();
    for n in [19, 23, 31, 37, 53, 83, 131] {
        module.push(json!({"n":n,"cyclotomic_coset_sizes":cyclotomic_cosets(n).map(|c|c.iter().map(Vec::len).collect::<Vec<_>>()),
            "stable_subspace_dimensions":stable_subspace_dimensions(n)}));
    }
    // A failure in one prerequisite must not discard earlier or unrelated
    // evidence. Preserve the complete stage outcome instead of returning early.
    let stage = |result: Result<Value, String>| match result {
        Ok(data) => json!({"status":"completed","data":data}),
        Err(error) => json!({"status":"error","error":error}),
    };
    let out = json!({"schema":"ecbench-ten-questions-controls/v2",
        "protocol":"CONTROLS-AMENDMENT-v2.md",
        "quasi_subfield":stage(quasi()),"nonlinear_bases":stage(nonlinear()),"degree_profiles":stage(degree_profiles()),
        "qubo":stage(qubo_control()),"covers":stage(covers()),"isogeny":stage(isogeny_controls()),
        "frobenius_module":module,"new_mathematical_theorem_claimed":false});
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_lib::cryptanalysis::pq_groebner_f2::F2BoolMono;
    #[test]
    fn qubo_distinguishes_constraint_violation_and_rejects_wrong_auxiliary() {
        let eq = F2BoolPoly::from_monos(
            vec![
                F2BoolMono { mask: 7 },
                F2BoolMono { mask: 1 },
                F2BoolMono { mask: 0 },
            ],
            3,
        );
        let q = Qubo::encode(std::slice::from_ref(&eq), 3);
        for a in 0..8 {
            let good = q.faithful_completion(a, 3);
            assert_eq!(q.energy(&good), i64::from(eval(&eq, a)));
        }
        assert!(and_table());
        // Exhaust every auxiliary choice in this small independent fixture.
        for a in 0..8 {
            let mut min = i64::MAX;
            for aux in 0..1usize << (q.vars - 3) {
                let bits: Vec<_> = (0..q.vars)
                    .map(|i| {
                        if i < 3 {
                            a >> i & 1 == 1
                        } else {
                            aux >> (i - 3) & 1 == 1
                        }
                    })
                    .collect();
                min = min.min(q.energy(&bits));
            }
            assert_eq!(min == 0, !eval(&eq, a));
        }
    }
    #[test]
    fn tiny_semaev_encoding_replays_exactly() {
        let out = qubo_control().unwrap();
        assert_eq!(
            out["faithful_completion_matches_violated_equation_count"],
            true
        );
        assert_eq!(out["original_assignments_examined"], 8192);
        assert_eq!(out["hardware_annealing"], false);
    }
}
