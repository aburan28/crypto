use crypto_lib::cryptanalysis::identity_certificate::{
    check_with_extra, issue, Expr, IdentityStatement, Verdict,
};
use std::fs::OpenOptions;
use std::io::Write;
fn main() {
    let output = std::env::args().nth(1).expect("new certificate output");
    let v = Expr::var;
    let plus = Expr::plus;
    let mul = Expr::times;
    let cert = issue(
        IdentityStatement {
            name: "two-limb exact product expansion".into(),
            source: "research/n83_speed_20261012/audit/ARITHMETIC_PROOF.md Proposition 1; elementary distributivity".into(),
            variables: ["a0", "a1", "b0", "b1", "w"].into_iter().map(str::to_owned).collect(),
            degree_bound: 4,
            lhs: mul(plus(v("a0"),mul(v("w"),v("a1"))),plus(v("b0"),mul(v("w"),v("b1")))),
            rhs: Expr::sum(vec![mul(v("a0"),v("b0")),mul(v("w"),plus(mul(v("a0"),v("b1")),mul(v("a1"),v("b0")))),Expr::product(vec![v("w"),v("w"),v("a1"),v("b1")])]),
        },
        202610126282,
        64,
        "Issued by independent speed-audit task from elementary distributivity, not copied from worker implementation; replayed by repository native certificate checker with 64 additional seed-separated points. This evaluated finite-field certificate does not prove integer overflow, Barrett quotient bounds, or end-to-end solver correctness.",
        None,
    ).expect("limb expansion holds");
    let verdict = check_with_extra(&cert, Some((202610126283,64))).expect("certificate well-formed");
    assert!(matches!(verdict,Verdict::Accept{..}),"{verdict:?}");
    let bytes = format!("{}\n",serde_json::to_string_pretty(&cert).unwrap());
    let mut f=OpenOptions::new().write(true).create_new(true).open(output).unwrap();
    f.write_all(bytes.as_bytes()).unwrap();
    f.sync_all().unwrap();
    println!("{{\"schema\":\"n83-limb-identity-audit/v1\",\"id\":\"{}\",\"status\":\"PASS\",\"certificate_points\":64,\"checker_extra_points\":64,\"checker_seed\":202610126283,\"formal_proof_pointer\":null,\"scope\":\"finite-field polynomial expansion only\"}}",cert.id);
}
