//! Optimal traversal strategies for isogeny chains (De Feo–Jao–Plût).
//! A strategy is a public, parameter-only object.

/// Optimal strategy for an n-step chain where one "multiplication" step
/// (point tripling / divisor doubling) costs `p` and one isogeny
/// evaluation costs `q`.  Returns the SIKE-style list of n-1 splits.
pub fn optimal(n: usize, p: f64, q: f64) -> Vec<usize> {
    if n <= 1 {
        return vec![];
    }
    let mut s: Vec<Vec<usize>> = vec![vec![]; n + 1];
    let mut c = vec![0f64; n + 1];
    for i in 2..=n {
        let mut best = (0usize, f64::INFINITY);
        for b in 1..i {
            let cost = c[i - b] + c[b] + (b as f64) * p + ((i - b) as f64) * q;
            if cost < best.1 {
                best = (b, cost);
            }
        }
        let b = best.0;
        let mut v = vec![b];
        v.extend(&s[i - b]);
        v.extend(&s[b]);
        s[i] = v;
        c[i] = best.1;
    }
    s[n].clone()
}

/// Generic strategy traversal.  `mult(x, s)` multiplies a kernel element by
/// l^s; `step(x)` builds an isogeny from an element of order l and returns
/// it; `eval(phi, x)` evaluates it.  `extra` are points pushed through every
/// step.  The control flow depends only on `n` and `strategy`.
pub fn traverse<K: Clone, I>(
    n: usize,
    strategy: &[usize],
    gen: K,
    extra: &mut [K],
    mut mult: impl FnMut(&K, usize) -> K,
    mut step: impl FnMut(&K) -> I,
    mut eval: impl FnMut(&I, &K) -> K,
) -> Vec<I> {
    assert!(n == 0 || strategy.len() == n - 1);
    let mut stack: Vec<(K, usize)> = vec![(gen, 0)];
    let mut idx = 0usize;
    let mut out = Vec::with_capacity(n);
    for row in 0..n {
        loop {
            let (q, k) = stack.last().unwrap().clone();
            if k >= n - 1 - row {
                break;
            }
            let s = strategy[idx];
            idx += 1;
            stack.push((mult(&q, s), k + s));
        }
        let (leaf, _) = stack.pop().unwrap();
        let phi = step(&leaf);
        for e in stack.iter_mut() {
            e.0 = eval(&phi, &e.0);
        }
        for x in extra.iter_mut() {
            *x = eval(&phi, x);
        }
        out.push(phi);
    }
    out
}
