use crypto_lib::{
    cryptanalysis::ghs_screen::{screen_curve, GhsCurveInput},
    hash::sha256,
};
use num_bigint::BigUint;
use num_traits::{One, Zero};
use serde_json::Value;
use std::{fs, path::PathBuf};

const EVIDENCE_DIR: &str = "research/sect113r1-weak-curve-20261006/addenda/ghs-three-node-20261007";
const PRIOR_CERTIFICATE: &str =
    "research/sect113r1-weak-curve-20261006/diagnostics/pre-admission-certificate.json";
const PRIOR_RAW_SHA256: &str = "8f77c0a21d5a51df372ce0a3a69784c58c3e404e5f1f52ad953c5e5714eeed24";
const PRIOR_SEMANTIC_SHA256: &str =
    "1cf4a3733e98c6acdc2219c24f29c15771fc67578bdb33dd426ab643e12572d8";
const PRODUCER_COMMIT: &str = "1f0d0bd7863c78ac9d7291b84e78ed937e8c4e2f";
const MODULUS: &str = "20000000000000000000000000201";
const CURVE_A: &str = "3088250ca6e7c7fe649ce85820f7";
const EXPECTED_GENUS: &str = "5192296858534827628530496329220095";
const EXPECTED_COVER_DEGREE: &str = "10384593717069655257060992658440192";

struct Case {
    slug: &'static str,
    role: &'static str,
    prior_curve_ref: &'static str,
    icv1: &'static str,
    curve_id: &'static str,
    curve_uid: &'static str,
    file: &'static str,
    b: &'static str,
    sha256: &'static str,
}

const CASES: &[Case] = &[
    Case {
        slug: "icv1-f2m113-tm122610772499221213-97df4ac6",
        role: "source (standard name: sect113r1)",
        prior_curve_ref: "sect113r1/source",
        icv1: "ICV1:f2m-113-99967757:-122610772499221213:10384593717069655379671765157661406:0x6942e38fc45c62366c09aa8204cd:unk:unk:r:97df4ac684cb",
        curve_id: "EC1N113Csect113r1hf529f17bd191",
        curve_uid: "urn:ec-record:1:sha256:f529f17bd1913792333a661e3557ad6b8e0ca2d4d02b939bc069d17d9fd94d97",
        file: "source.json",
        b: "e8bee4d3e2260744188be0e9c723",
        sha256: "08c763d46c87c8c416fd6688d99d9b84a40c2a256ecf31127e06a8eefb458b5d",
    },
    Case {
        slug: "icv1-f2m113-tm122610772499221213-fd54d0eb",
        role: "degree-5 codomain A",
        prior_curve_ref: "sect113r1/degree5-a",
        icv1: "ICV1:f2m-113-99967757:-122610772499221213:10384593717069655379671765157661406:0xb1df3d9c423aa919217735a6d6ba:unk:unk:r:fd54d0ebcd45",
        curve_id: "EC1N113Crbh921ab2cd913f",
        curve_uid: "urn:ec-record:1:sha256:921ab2cd913f1d63106c1a55555cd869ae110bff6992dd23ee2522a4b97616c6",
        file: "degree5-a.json",
        b: "109267245489e254e8f14002629a1",
        sha256: "4bc30b46310dffd334ba858bbe6d67cec394c556a69e7e85a2a7e3fb00fad88f",
    },
    Case {
        slug: "icv1-f2m113-tm122610772499221213-5de3030f",
        role: "degree-5 codomain B",
        prior_curve_ref: "sect113r1/degree5-b",
        icv1: "ICV1:f2m-113-99967757:-122610772499221213:10384593717069655379671765157661406:0x17ef6a3b098891f6b4c529b62f0dd:unk:unk:r:5de3030fefcc",
        curve_id: "EC1N113Crbh2ee3581888f2",
        curve_uid: "urn:ec-record:1:sha256:2ee3581888f2dcbc20ded6386aac140f0ced5d5b8ff92aab74887ecd19f19386",
        file: "degree5-b.json",
        b: "162b1a595685a1387c82647bf44cf",
        sha256: "79fd42d92e984cc4de339fe245b2d3aa2f247aa8536777975570bf241919be6a",
    },
];

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn parse_hex(value: &str) -> BigUint {
    BigUint::parse_bytes(value.as_bytes(), 16).expect("valid frozen hexadecimal integer")
}

fn sha256_hex(bytes: &[u8]) -> String {
    hex::encode(sha256(bytes))
}

fn pell_power(exponent: u32) -> (BigUint, BigUint) {
    // Exact recurrence for (3 + 2*sqrt(2))^exponent = A + B*sqrt(2).
    let mut a = BigUint::one();
    let mut b = BigUint::zero();
    for _ in 0..exponent {
        let next_a = &a * 3u8 + &b * 4u8;
        let next_b = &a * 2u8 + &b * 3u8;
        a = next_a;
        b = next_b;
    }
    (a, b)
}

fn less_than_a_plus_b_sqrt_two(value: &BigUint, a: &BigUint, b: &BigUint) -> bool {
    if value <= a {
        return true;
    }
    let difference = value - a;
    difference.pow(2) < b.pow(2) * 2u8
}

#[test]
fn preserved_outputs_match_hashes_and_recompute_exactly() {
    let evidence = root().join(EVIDENCE_DIR);
    let manifest_bytes = fs::read(evidence.join("manifest.json")).expect("read manifest");
    let manifest: Value = serde_json::from_slice(&manifest_bytes).expect("parse manifest");
    assert_eq!(
        manifest["schema"],
        "crypto.sect113r1-ghs-three-node-screen/v1"
    );
    assert_eq!(
        manifest["evidence_status"],
        "ADDITIVE_PRE_ADMISSION_STRUCTURAL_DIAGNOSTIC"
    );
    assert_eq!(manifest["admitted_scientific_run"], false);
    assert_eq!(manifest["independent_implementation_validation"], false);
    assert_eq!(manifest["producer"]["git_commit"], PRODUCER_COMMIT);
    assert_eq!(manifest["producer"]["profile"], "release");
    assert_eq!(
        manifest["producer"]["cargo"],
        "cargo 1.99.0 (5f94df478 2026-08-27)"
    );
    assert_eq!(
        manifest["producer"]["rustc"],
        "rustc 1.99.0 (b940084d7 2026-09-28)"
    );
    assert_eq!(
        manifest["producer"]["binary"]["sha256"],
        "498b711179ec8a9dd9e7a12528b2025b1b7f7a89f210c28ccc6f625dc899f9ca"
    );
    assert_eq!(manifest["producer"]["binary"]["bytes"], 2_402_944);

    let predecessor = &manifest["frozen_predecessor"];
    assert_eq!(
        predecessor["path"],
        "../../diagnostics/pre-admission-certificate.json"
    );
    assert_eq!(predecessor["raw_sha256"], PRIOR_RAW_SHA256);
    assert_eq!(predecessor["semantic_sha256"], PRIOR_SEMANTIC_SHA256);
    assert_eq!(predecessor["modified_by_this_addendum"], false);
    assert_eq!(
        predecessor["prior_ghs_status"],
        "NOT_CERTIFIED_MAGIC_NUMBER_UNCOMPUTED"
    );
    let prior_from_manifest = evidence.join(
        predecessor["path"]
            .as_str()
            .expect("predecessor path is a string"),
    );
    assert_eq!(
        fs::canonicalize(&prior_from_manifest).expect("canonical predecessor path"),
        fs::canonicalize(root().join(PRIOR_CERTIFICATE)).expect("canonical frozen path")
    );
    let prior_bytes = fs::read(prior_from_manifest).expect("read prior certificate");
    assert_eq!(sha256_hex(&prior_bytes), PRIOR_RAW_SHA256);
    let prior: Value = serde_json::from_slice(&prior_bytes).expect("parse prior certificate");
    assert_eq!(prior["certificate_sha256"], PRIOR_SEMANTIC_SHA256);
    assert_eq!(
        prior["source"]["ghs_screen_status"],
        "NOT_CERTIFIED_MAGIC_NUMBER_UNCOMPUTED"
    );

    for source in manifest["producer"]["source_files"]
        .as_array()
        .expect("source-files array")
    {
        let path = source["path"].as_str().expect("source path string");
        let expected = source["sha256"].as_str().expect("source hash string");
        let bytes = fs::read(root().join(path)).expect("read producer source");
        assert_eq!(sha256_hex(&bytes), expected, "producer source {path}");
    }

    assert_eq!(manifest["protocol"]["absolute_degree"], 113);
    assert_eq!(
        manifest["protocol"]["field_modulus"],
        format!("0x{MODULUS}")
    );
    assert_eq!(manifest["protocol"]["a"], format!("0x{CURVE_A}"));
    assert_eq!(manifest["protocol"]["genus_bound"], "64");

    let expected_genus = BigUint::one() << 112usize;
    let expected_genus = expected_genus - BigUint::one();
    let expected_cover_degree = BigUint::one() << 113usize;
    let executions = manifest["executions"]
        .as_array()
        .expect("manifest executions array");
    assert_eq!(executions.len(), CASES.len());

    for case in CASES {
        let bytes = fs::read(evidence.join(case.file)).expect("read preserved output");
        assert_eq!(sha256_hex(&bytes), case.sha256, "{} hash", case.slug);
        let json: Value = serde_json::from_slice(&bytes).expect("parse preserved output");
        assert_eq!(json["schema"], "ghs-screen/v1");
        assert_eq!(json["curve"]["model"], "y^2+xy=x^3+a*x^2+b");
        assert_eq!(json["curve"]["absolute_degree"], 113);
        assert_eq!(json["curve"]["modulus"], format!("0x{MODULUS}"));
        assert_eq!(json["curve"]["a"], format!("0x{CURVE_A}"));
        assert_eq!(json["curve"]["b"], format!("0x{}", case.b));
        assert_eq!(json["summary"]["factorisations"], 1);
        assert_eq!(json["summary"]["candidates_within_bound"], 0);
        assert_eq!(json["summary"]["genus_bound"], "64");
        assert_eq!(json["summary"]["best_relative_degree"], 113);
        assert_eq!(json["summary"]["best_base_degree"], 1);
        assert_eq!(json["summary"]["best_magic_number"], 113);
        assert_eq!(json["summary"]["best_genus"], EXPECTED_GENUS);
        assert_eq!(json["towers"].as_array().expect("tower array").len(), 1);
        assert_eq!(json["towers"][0]["relative_degree"], 113);
        assert_eq!(json["towers"][0]["base_degree"], 1);
        assert_eq!(
            json["towers"][0]["field_tower"],
            "F_2 <= F_(2^1) <= F_(2^113)"
        );
        assert_eq!(json["towers"][0]["magic_number"], 113);
        assert_eq!(json["towers"][0]["genus"], EXPECTED_GENUS);
        assert_eq!(json["towers"][0]["ghs_type"], "I");
        assert_eq!(
            json["towers"][0]["artin_schreier_cover_degree"],
            EXPECTED_COVER_DEGREE
        );
        assert_eq!(json["towers"][0]["within_genus_bound"], false);

        let input = GhsCurveInput {
            absolute_degree: 113,
            modulus: parse_hex(MODULUS),
            a: parse_hex(CURVE_A),
            b: parse_hex(case.b),
        };
        let report = screen_curve(&input, &BigUint::from(64u8)).expect("exact GHS screen");
        assert_eq!(report.rows.len(), 1, "{} factorization count", case.slug);
        let row = &report.rows[0];
        assert_eq!(row.relative_degree, 113);
        assert_eq!(row.base_degree, 1);
        assert_eq!(row.magic_number, 113);
        assert!(row.type_i);
        assert_eq!(row.genus, expected_genus);
        assert_eq!(row.cover_degree, expected_cover_degree);
        assert!(!row.within_genus_bound);

        let execution = executions
            .iter()
            .find(|execution| execution["node"] == case.slug)
            .expect("manifest execution for node");
        assert_eq!(execution["role"], case.role);
        assert_eq!(execution["prior_curve_ref"], case.prior_curve_ref);
        assert_eq!(execution["icv1"], case.icv1);
        assert_eq!(execution["icv1_slug"], case.slug);
        assert_eq!(execution["curve_id"], case.curve_id);
        assert_eq!(execution["curve_uid"], case.curve_uid);
        assert_eq!(execution["b"], format!("0x{}", case.b));
        assert_eq!(execution["output"], case.file);
        assert_eq!(execution["output_sha256"], case.sha256);
        assert_eq!(execution["replay_output_sha256"], case.sha256);
        assert_eq!(execution["replay_byte_equal"], true);
        assert_eq!(execution["first_exit_code"], 0);
        assert_eq!(execution["replay_exit_code"], 0);
        assert_eq!(execution["result"]["factorizations"], 1);
        assert_eq!(execution["result"]["relative_degree"], 113);
        assert_eq!(execution["result"]["base_degree"], 1);
        assert_eq!(execution["result"]["magic_number"], 113);
        assert_eq!(execution["result"]["ghs_type"], "I");
        assert_eq!(execution["result"]["genus"], EXPECTED_GENUS);
        assert_eq!(
            execution["result"]["artin_schreier_cover_degree"],
            EXPECTED_COVER_DEGREE
        );
        assert_eq!(execution["result"]["within_genus_bound"], false);
    }

    let mut expected_csv =
        "icv1_slug,magic_number,ghs_type,genus,genus_bound,within_bound\n".to_owned();
    for case in CASES {
        expected_csv.push_str(&format!("{},113,I,{EXPECTED_GENUS},64,false\n", case.slug));
    }
    assert_eq!(
        fs::read_to_string(evidence.join("ghs-magic-comparison.csv"))
            .expect("read quantitative visual data"),
        expected_csv
    );
    for visual in ["ghs-magic-comparison.svg", "ghs-screen-flow.dot"] {
        let text = fs::read_to_string(evidence.join(visual)).expect("read visual source");
        for case in CASES {
            assert!(text.contains(case.slug), "{visual} omits {}", case.slug);
        }
    }
}

#[test]
fn weil_capacity_gate_excludes_genus_at_most_44() {
    let subgroup_order = BigUint::parse_bytes(b"5192296858534827689835882578830703", 10)
        .expect("valid subgroup order");

    // #Jac(C)(F_2) <= (1 + sqrt(2))^(2g) = (3 + 2*sqrt(2))^g.
    // Compare A+B*sqrt(2) with n by squaring the positive residual.
    let (a44, b44) = pell_power(44);
    assert!(!less_than_a_plus_b_sqrt_two(&subgroup_order, &a44, &b44));

    // The capacity bound crosses n between 44 and 45. This does not prove
    // that a genus-45 Jacobian has the subgroup; it only removes the size
    // obstruction at that point.
    let (a45, b45) = pell_power(45);
    assert!(less_than_a_plus_b_sqrt_two(&subgroup_order, &a45, &b45));
}
