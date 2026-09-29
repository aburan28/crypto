mod profiled {
    use super::*;
    include!("profiled_quadratic.rs");
}
fn main() {
    for n in [16, 24] {
        for seed in [17, 937] {
            for family in ["planted", "cross_planted", "unplanted"] {
                let (system, _) = fixture(n, seed, family);
                let expected = solve(&system, n, "quadratic_state", 200000);
                for rep in 0..3 {
                    let tick = Instant::now();
                    let (got, p) = profiled::measure(&system, n, 200000);
                    let total = tick.elapsed().as_nanos();
                    assert_eq!(got.outcome, expected.outcome);
                    assert_eq!(got.logical, expected.logical);
                    assert_eq!(got.trace, expected.trace);
                    if let Outcome::Sat(model) = got.outcome {
                        assert!(satisfies(&system, model));
                    }
                    let measured =
                        p.compile_ns + p.trace_ns + p.affine_ns + p.decision_ns + p.specialize_ns;
                    assert!(measured <= total);
                    println!("{{\"n\":{n},\"seed\":{seed},\"family\":\"{family}\",\"rep\":{rep},\"nodes\":{},\"trace\":{},\"compile_ns\":{},\"trace_ns\":{},\"affine_ns\":{},\"decision_ns\":{},\"specialize_ns\":{},\"total_ns\":{total},\"remainder_ns\":{},\"matched\":true}}",got.logical.nodes,got.trace,p.compile_ns,p.trace_ns,p.affine_ns,p.decision_ns,p.specialize_ns,total-measured);
                }
            }
        }
    }
}
