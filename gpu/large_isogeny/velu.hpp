#pragma once

// Bounded prime-field prototype. p < 2^31 so products fit uint64_t.
#include <cstdint>

#ifdef __CUDACC__
#define ISO_HD __host__ __device__
#else
#define ISO_HD
#endif

namespace iso {
struct Curve { uint64_t p, a, b; };
struct Point { uint64_t x, y; bool infinity; };

ISO_HD inline uint64_t add(uint64_t a, uint64_t b, uint64_t p) { return (a + b) % p; }
ISO_HD inline uint64_t sub(uint64_t a, uint64_t b, uint64_t p) { return (a + p - b) % p; }
ISO_HD inline uint64_t mul(uint64_t a, uint64_t b, uint64_t p) { return (a * b) % p; }
ISO_HD inline uint64_t power(uint64_t a, uint64_t n, uint64_t p) {
    uint64_t r = 1;
    for (; n; n >>= 1, a = mul(a, a, p)) if (n & 1) r = mul(r, a, p);
    return r;
}
// p is prime and a != 0. The caller checks the exceptional denominator.
ISO_HD inline uint64_t inverse(uint64_t a, uint64_t p) { return power(a, p - 2, p); }
ISO_HD inline Point infinity() { return {0, 0, true}; }
ISO_HD inline bool equal(Point a, Point b) {
    return a.infinity == b.infinity && (a.infinity || (a.x == b.x && a.y == b.y));
}
ISO_HD inline Point negate(Point a, Curve e) {
    if (!a.infinity) a.y = sub(0, a.y, e.p);
    return a;
}
ISO_HD inline bool on_curve(Point q, Curve e) {
    if (q.infinity) return true;
    const auto lhs = mul(q.y, q.y, e.p);
    const auto rhs = add(add(mul(mul(q.x, q.x, e.p), q.x, e.p), mul(e.a, q.x, e.p), e.p), e.b, e.p);
    return lhs == rhs;
}
ISO_HD inline Point sum(Point x, Point y, Curve e) {
    if (x.infinity) return y;
    if (y.infinity) return x;
    uint64_t num, den;
    if (x.x == y.x) {
        if (x.y != y.y || x.y == 0) return infinity();
        num = add(mul(3, mul(x.x, x.x, e.p), e.p), e.a, e.p);
        den = mul(2, x.y, e.p);
    } else {
        num = sub(y.y, x.y, e.p);
        den = sub(y.x, x.x, e.p);
    }
    const auto slope = mul(num, inverse(den, e.p), e.p);
    const auto rx = sub(sub(mul(slope, slope, e.p), x.x, e.p), y.x, e.p);
    return {rx, sub(mul(slope, sub(x.x, rx, e.p), e.p), x.y, e.p), false};
}
ISO_HD inline Point scalar(Point p, uint64_t n, Curve e) {
    Point r = infinity();
    for (; n; n >>= 1, p = sum(p, p, e)) if (n & 1) r = sum(r, p, e);
    return r;
}
// Full nonzero kernel in pairs Q,-Q. Caller certifies distinctness/order.
ISO_HD inline Point velu(Point p, const Point *kernel, uint64_t count, Curve e) {
    if (p.infinity) return p;
    uint64_t x = p.x, y = p.y;
    for (uint64_t i = 0; i < count; ++i) {
        const Point s = sum(p, kernel[i], e);
        if (s.infinity) return infinity(); // p belongs to the kernel
        x = add(x, sub(s.x, kernel[i].x, e.p), e.p);
        y = add(y, sub(s.y, kernel[i].y, e.p), e.p);
    }
    return {x, y, false};
}
ISO_HD inline Curve codomain(const Point *kernel, uint64_t count, Curve e) {
    uint64_t v = 0, w = 0;
    for (uint64_t i = 0; i < count; ++i) {
        const auto x = kernel[i].x;
        v = add(v, add(mul(3, mul(x, x, e.p), e.p), e.a, e.p), e.p);
        w = add(w, add(add(mul(5, mul(mul(x, x, e.p), x, e.p), e.p),
                            mul(3, mul(e.a, x, e.p), e.p), e.p), mul(2, e.b, e.p), e.p), e.p);
    }
    return {e.p, sub(e.a, mul(5, v, e.p), e.p), sub(e.b, mul(7, w, e.p), e.p)};
}
} // namespace iso
