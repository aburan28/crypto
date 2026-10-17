//! Exact value-level local transitions, bounded beam proposals. The DDT product
//! is a random-round-key model score, not a fixed-key probability or hull bound.
use crate::cipher::{linear, Block, Config, SBOX};
use crate::{diagnostics::ddt, hex};
use serde::Serialize;

#[derive(Clone, Debug, Serialize)]
pub struct RoundStep {
    pub round: usize,
    pub before_sub_bytes: String,
    pub after_sub_bytes: String,
    pub after_linear: String,
    pub active_sboxes: usize,
    pub ddt_weight: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct Trail {
    pub output_difference: String,
    pub ddt_product_weight: f64,
    pub steps: Vec<RoundStep>,
}

#[derive(Serialize)]
pub struct TrailReport {
    pub schema: String,
    pub cipher: Config,
    pub input_difference: String,
    pub beam_width: usize,
    pub branch_cap_per_sbox: usize,
    pub expanded_local_transitions: u64,
    pub discarded_partial_nodes: u64,
    pub trails: Vec<Trail>,
    pub scope: String,
}

#[derive(Clone)]
struct Node {
    state: Block,
    weight: f64,
    steps: Vec<RoundStep>,
}

pub fn search(
    config: Config,
    input: Block,
    beam: usize,
    branch: usize,
) -> Result<TrailReport, String> {
    config.check()?;
    if input == [0; 16] {
        return Err("zero input difference is excluded".into());
    }
    if !(1..=256).contains(&beam) || !(1..=127).contains(&branch) {
        return Err("beam must be 1..256; branch cap 1..127".into());
    }
    let table = ddt(&SBOX);
    let rows: Vec<Vec<(u8, f64)>> = (0..256)
        .map(|a| {
            let mut v: Vec<_> = (0..256)
                .filter(|&b| table[a][b] > 0)
                .map(|b| (b as u8, -(table[a][b] as f64 / 256.0).log2()))
                .collect();
            v.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
            v
        })
        .collect();
    let mut frontier = vec![Node {
        state: input,
        weight: 0.0,
        steps: vec![],
    }];
    let mut expanded = 0;
    let mut discarded = 0;
    for round in 1..=config.rounds {
        let mut next = Vec::new();
        for node in frontier {
            let mut partial = vec![(node.state, 0.0)];
            for pos in 0..16 {
                let mut out = Vec::new();
                let row = &rows[node.state[pos] as usize];
                for (s, w) in partial {
                    for &(b, cost) in row.iter().take(branch) {
                        let mut t = s;
                        t[pos] = b;
                        out.push((t, w + cost));
                        expanded += 1;
                    }
                }
                out.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
                discarded += out.len().saturating_sub(beam) as u64;
                out.truncate(beam);
                partial = out;
            }
            for (after_sub, w) in partial {
                let mut state = after_sub;
                linear(
                    &mut state,
                    round < config.rounds || config.final_mix_columns,
                );
                let mut steps = node.steps.clone();
                steps.push(RoundStep {
                    round,
                    before_sub_bytes: hex(&node.state),
                    after_sub_bytes: hex(&after_sub),
                    after_linear: hex(&state),
                    active_sboxes: node.state.iter().filter(|&&x| x != 0).count(),
                    ddt_weight: w,
                });
                next.push(Node {
                    state,
                    weight: node.weight + w,
                    steps,
                });
            }
        }
        next.sort_by(|a, b| a.weight.total_cmp(&b.weight).then(a.state.cmp(&b.state)));
        discarded += next.len().saturating_sub(beam) as u64;
        next.truncate(beam);
        frontier = next;
    }
    Ok(TrailReport {schema:"aes-discover/trails/v1".into(),cipher:config,input_difference:hex(&input),
        beam_width:beam,branch_cap_per_sbox:branch,expanded_local_transitions:expanded,
        discarded_partial_nodes:discarded,
        trails:frontier.into_iter().map(|n|Trail {output_difference:hex(&n.state),ddt_product_weight:n.weight,steps:n.steps}).collect(),
        scope:"heuristic value-level trail proposals only; local DDT transitions are exact; no feasibility, fixed-key probability, differential-hull probability, optimality or absence certificate".into()})
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::block;
    #[test]
    fn one_active_byte_forces_one_four_sixteen() {
        let r = search(
            Config {
                key_bits: 128,
                rounds: 3,
                final_mix_columns: true,
            },
            block("01000000000000000000000000000000").unwrap(),
            4,
            2,
        )
        .unwrap();
        for t in r.trails {
            assert_eq!(
                t.steps.iter().map(|s| s.active_sboxes).collect::<Vec<_>>(),
                [1, 4, 16]
            );
            assert_eq!(t.steps[0].ddt_weight, 6.0);
        }
    }
}
