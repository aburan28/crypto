// Sanity-check binary: runs Buchberger with each built-in strategy on cyclic-n
// and prints step / basis / arith-op counts. This is the baseline table any
// learned policy has to beat.
//
//     cargo run --release --bin demo

use gbrl::env::Env;
use gbrl::semaev::cyclic;
use gbrl::strategy::{NormalStrategy, RandomStrategy, Strategy, SugarStrategy};
use rand::SeedableRng;
use rand::rngs::StdRng;

fn run(strategy: &mut dyn Strategy, n: usize) {
    let polys = cyclic(n);
    let mut env = Env::new(polys, 200_000);
    while !env.done() {
        let idx = strategy.select(&env.state);
        env.step(idx);
    }
    let obs = env.observe();
    println!(
        "  {:>7}  cyclic-{}  steps={:<6} basis={:<4} arith_ops={}",
        strategy.name(), n, obs.step_count, obs.basis_size, obs.arith_ops
    );
}

fn main() {
    for n in 3..=5 {
        run(&mut SugarStrategy, n);
        run(&mut NormalStrategy, n);
        run(&mut RandomStrategy(StdRng::seed_from_u64(0xC0FFEE)), n);
        println!();
    }
}
