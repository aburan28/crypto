# Generates x86-64 (Intel syntax) Montgomery CIOS multiplication with MULX/ADCX/ADOX for N limbs.
# Buffer layout at [rsi]: a[N] | b[N] | p[N] | pinv | out[N+1]
import sys
def gen(N):
    A, B, P = 0, 8*N, 16*N
    PINV = 24*N
    OUT = 24*N + 8
    if N == 8:
        pool = ["rbp", "rcx", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15"]
        lo, hi, Z = "rax", "rbx", "rdi"
        save = ["rbx", "rbp"]
    elif N == 4:
        pool = ["rcx", "r8", "r9", "r10", "r11", "r12"]
        lo, hi, Z = "rax", "r13", "rdi"
        save = []
    else:
        raise SystemExit("N")
    L = []
    for r in save:
        L.append(f"push {r}")
    R = pool[:N+1]
    X = pool[N+1]
    for r in R:
        L.append(f"xor {r:s}, {r:s}")
    for i in range(N):
        # t += a * b_i
        L.append(f"mov rdx, qword ptr [rsi + {B + 8*i}]")
        L.append(f"xor {Z}, {Z}")
        for j in range(N):
            L.append(f"mulx {hi}, {lo}, qword ptr [rsi + {A + 8*j}]")
            L.append(f"adox {R[j]}, {lo}")
            L.append(f"adcx {R[j+1]}, {hi}")
        L.append(f"mov {X}, 0")
        L.append(f"adox {R[N]}, {Z}")
        L.append(f"adcx {X}, {Z}")
        L.append(f"adox {X}, {Z}")
        # m = t0 * pinv; t += m * p; shift
        L.append(f"mov rdx, {R[0]}")
        L.append(f"imul rdx, qword ptr [rsi + {PINV}]")
        L.append(f"xor {Z}, {Z}")
        for j in range(N):
            L.append(f"mulx {hi}, {lo}, qword ptr [rsi + {P + 8*j}]")
            L.append(f"adox {R[j]}, {lo}")
            L.append(f"adcx {R[j+1]}, {hi}")
        L.append(f"adox {R[N]}, {Z}")
        L.append(f"adcx {X}, {Z}")
        L.append(f"adox {X}, {Z}")
        # rotate: new t_j = R[j+1], new t_N = X, new X = old R[0] (zero)
        R, X = R[1:] + [X], R[0]
    for j in range(N+1):
        L.append(f"mov qword ptr [rsi + {OUT + 8*j}], {R[j]}")
    for r in reversed(save):
        L.append(f"pop {r}")
    return L
for N in (4, 8):
    lines = gen(N)
    print(f"// N = {N}: {len(lines)} instructions")
    print(f"pub const MONT_ADX_{N}: &str = concat!(")
    for l in lines:
        print(f'    "{l}\\n",')
    print(");")
