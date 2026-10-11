//! Generator of the x86-64 MULX/ADCX/ADOX Montgomery multiplication (CIOS, accumulator in
//! registers) embedded in `fpm::adx`, for N = 4 and N = 8 limbs. Intel syntax; buffer at [rsi]:
//! a[N] | b[N] | p[N] | pinv | out[N+1]. The native replacement of the earlier
//! `scripts/gen_mont_adx.py`, emitting the same text (`tests/adx_gen.rs`).

/// The instructions for N limbs (N = 4 or 8).
pub fn mont_adx(n: usize) -> Vec<String> {
    let (a, b, p) = (0, 8 * n, 16 * n);
    let pinv = 24 * n;
    let out = 24 * n + 8;
    let (pool, lo, hi, z, save): (Vec<&str>, &str, &str, &str, Vec<&str>) = match n {
        8 => (
            vec![
                "rbp", "rcx", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15",
            ],
            "rax",
            "rbx",
            "rdi",
            vec!["rbx", "rbp"],
        ),
        4 => (
            vec!["rcx", "r8", "r9", "r10", "r11", "r12"],
            "rax",
            "r13",
            "rdi",
            vec![],
        ),
        _ => panic!("mont_adx: N must be 4 or 8"),
    };
    let mut l: Vec<String> = save.iter().map(|r| format!("push {r}")).collect();
    let mut r: Vec<&str> = pool[..n + 1].to_vec();
    let mut x = pool[n + 1];
    for reg in &r {
        l.push(format!("xor {reg}, {reg}"));
    }
    for i in 0..n {
        // t += a * b_i
        l.push(format!("mov rdx, qword ptr [rsi + {}]", b + 8 * i));
        l.push(format!("xor {z}, {z}"));
        for j in 0..n {
            l.push(format!("mulx {hi}, {lo}, qword ptr [rsi + {}]", a + 8 * j));
            l.push(format!("adox {}, {lo}", r[j]));
            l.push(format!("adcx {}, {hi}", r[j + 1]));
        }
        l.push(format!("mov {x}, 0"));
        l.push(format!("adox {}, {z}", r[n]));
        l.push(format!("adcx {x}, {z}"));
        l.push(format!("adox {x}, {z}"));
        // m = t0 * pinv; t += m * p; shift
        l.push(format!("mov rdx, {}", r[0]));
        l.push(format!("imul rdx, qword ptr [rsi + {pinv}]"));
        l.push(format!("xor {z}, {z}"));
        for j in 0..n {
            l.push(format!("mulx {hi}, {lo}, qword ptr [rsi + {}]", p + 8 * j));
            l.push(format!("adox {}, {lo}", r[j]));
            l.push(format!("adcx {}, {hi}", r[j + 1]));
        }
        l.push(format!("adox {}, {z}", r[n]));
        l.push(format!("adcx {x}, {z}"));
        l.push(format!("adox {x}, {z}"));
        // rotate: new t_j = r[j+1], new t_N = x, new x = old r[0] (zero)
        let r0 = r.remove(0);
        r.push(x);
        x = r0;
    }
    for (j, reg) in r.iter().enumerate() {
        l.push(format!("mov qword ptr [rsi + {}], {reg}", out + 8 * j));
    }
    for reg in save.iter().rev() {
        l.push(format!("pop {reg}"));
    }
    l
}

/// The Rust source text for both routines, as `gen_mont_adx` prints it.
pub fn source() -> String {
    let mut s = String::new();
    for n in [4, 8] {
        let lines = mont_adx(n);
        s.push_str(&format!(
            "// N = {n}: {} instructions\npub const MONT_ADX_{n}: &str = concat!(\n",
            lines.len()
        ));
        for line in &lines {
            s.push_str(&format!("    \"{line}\\n\",\n"));
        }
        s.push_str(");\n");
    }
    s
}
