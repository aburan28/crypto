"""Exact per-target changes for a separately declared batch tournament."""


def lazy_descent(text):
    start = text.index('    fn solve_by_walking(')
    end = text.index('    pub fn solve(&self, q: &BinaryPoint)', start)
    body = text[start:end]
    needle = '        let stride_point = fc.mul_u64(*g, stride);'
    prefix = '''        if report.trials >= self.opts.max_trials { return Some(None); }
        let initial = fc.add(fc.mul_u64(*g, a0), fc.mul_u64(q_fast, b));
        // Many targets descend at the first state. Pay for that probe before
        // initializing the other 63 walks; keep the same coefficients/order.
        report.trials += 1;
        if initial.infinity {
            if let Some(d) = solve_for_d(&BigUint::from(a0), &BigUint::from(b), &self.kc.subgroup_order) {
                if self.kc.mul(self.kc.generator(), &d) == *q { return Some(Some(d)); }
            }
        } else if let Some(idxs) = pair.decompose_fast(initial, m) {
            if let Some(d) = self.logarithm_from(q, &idxs, a0, b) { return Some(Some(d)); }
        }
        if report.trials >= self.opts.max_trials { return Some(None); }
'''
    assert body.count(needle) == 1
    body = body.replace(needle, prefix + needle, 1)
    body = body.replace('let mut state = fc.add(fc.mul_u64(*g, a0), fc.mul_u64(q_fast, b));',
                        'let mut state = initial;', 1)
    body = body.replace('        while report.trials < self.opts.max_trials {',
                        '        let mut skip_first = true;\n        while report.trials < self.opts.max_trials {', 1)
    needle = '            for (state, (a, b)) in states.iter().zip(coefficients.iter()) {\n'
    assert body.count(needle) == 1
    body = body.replace(needle, needle + '''                if skip_first { skip_first = false; continue; }
                if report.trials >= self.opts.max_trials { return Some(None); }
''', 1)
    return text[:start] + body + text[end:]


def fast_descent_check(text):
    start = text.index("impl<'a> IndividualLogSolver<'a> {")
    end = text.index('/// Outcome of one descent probe.', start)
    body = text[start:end]
    body = body.replace('kc.mul(kc.generator(), &d) == *q', 'self.verify_log(q, &d)')
    # The replacement above also matches inside self.kc; fix that prefix.
    body = body.replace('self.self.verify_log(q, &d)', 'self.verify_log(q, &d)')
    body = body.replace('kc.mul(g, &d) == *q', 'self.verify_log(q, &d)')
    body = body.replace('`[d]G = Q` in the general arithmetic before it is returned.',
                        '`[d]G = Q` in the curve arithmetic before it is returned.')
    at = body.index('\n') + 1
    helper = '''    fn verify_log(&self, q: &BinaryPoint, d: &BigUint) -> bool {
        if d >= &self.kc.subgroup_order { return false; }
        if let Some((fc, g)) = &self.fast {
            fc.mul(*g, d) == fc.lift(q)
        } else {
            self.kc.mul(self.kc.generator(), d) == *q
        }
    }

'''
    body = body[:at] + helper + body[at:]
    return text[:start] + body + text[end:]
