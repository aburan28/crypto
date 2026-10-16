// Bit-for-bit oracle for walk_step_problem.py: one canonical walk step,
// P + sigma^j(P) with j = 3 + ((HW(x) >> 1) & 7), through include/ref.h's
// Ref<CfgF131> (ScalarOnb field) and the device's special-case-free addPtRaw.
// Reads "x0 x1 x2 y0 y1 y2" hex limbs per line on stdin, writes the six
// result limbs per line.  Built and driven by check_reference.py.
#include <cstdio>
#include "../include/curveparams.h"

typedef Ref<CfgF131> S;

int main() {
    unsigned long long in[6];
    while (scanf("%llx %llx %llx %llx %llx %llx", &in[0], &in[1], &in[2], &in[3], &in[4], &in[5]) == 6) {
        const S::Point p = S::make(S::fromLimbs(in), S::fromLimbs(in + 3));
        const S::Point q = S::addPtRaw(p, S::frob(p, S::jOf(S::weight(p.x))));
        printf("%llx %llx %llx %llx %llx %llx\n", q.x.v[0], q.x.v[1], q.x.v[2], q.y.v[0], q.y.v[1], q.y.v[2]);
    }
    return 0;
}
