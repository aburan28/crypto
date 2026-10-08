#include "velu.hpp"
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

using namespace iso;

static bool prime(uint64_t n) {
    if (n < 2) return false;
    for (uint64_t d = 2; d * d <= n; ++d) if (n % d == 0) return false;
    return true;
}

static uint64_t count_points(Curve e) {
    uint64_t n = 1;
    for (uint64_t x = 0; x < e.p; ++x) {
        const uint64_t rhs = add(add(mul(mul(x, x, e.p), x, e.p), mul(e.a, x, e.p), e.p), e.b, e.p);
        n += rhs == 0 ? 1 : (power(rhs, (e.p - 1) / 2, e.p) == 1 ? 2 : 0);
    }
    return n;
}

static Point find_generator(Curve e, uint64_t cofactor, uint64_t ell) {
    for (uint64_t x = 0; x < e.p; ++x) {
        for (uint64_t y = 0; y < e.p; ++y) {
            Point p{x, y, false};
            if (!on_curve(p, e)) continue;
            p = scalar(p, cofactor, e);
            if (!p.infinity && scalar(p, ell, e).infinity) return p;
        }
    }
    return infinity();
}

int main() {
    constexpr uint64_t p = 1009, ell = 101;
    static_assert(ell < p);
    assert(prime(p) && prime(ell));
    Curve e{};
    uint64_t order = 0;
    for (uint64_t a = 0; a < 30 && !order; ++a) {
        for (uint64_t b = 1; b < 30 && !order; ++b) {
            Curve candidate{p, a, b};
            if (add(mul(4, mul(mul(a, a, p), a, p), p), mul(27, mul(b, b, p), p), p) == 0) continue;
            uint64_t n = count_points(candidate);
            if (n % ell == 0) { e = candidate; order = n; }
        }
    }
    assert(order != 0);
    const Point generator = find_generator(e, order / ell, ell);
    assert(!generator.infinity);
    std::vector<Point> kernel;
    for (Point q = generator; !q.infinity; q = sum(q, generator, e)) kernel.push_back(q);
    assert(kernel.size() == ell - 1);
    for (size_t i = 0; i < kernel.size(); ++i) {
        assert(on_curve(kernel[i], e));
        for (size_t j = 0; j < i; ++j) assert(!equal(kernel[i], kernel[j]));
        assert(velu(kernel[i], kernel.data(), kernel.size(), e).infinity);
    }
    const Curve target = codomain(kernel.data(), kernel.size(), e);
    assert(count_points(target) == order);
    uint64_t checked = 0;
    for (uint64_t x = 0; x < p && checked < 30; ++x) {
        for (uint64_t y = 0; y < p && checked < 30; ++y) {
            Point q{x, y, false};
            if (!on_curve(q, e)) continue;
            const Point r = velu(q, kernel.data(), kernel.size(), e);
            assert(on_curve(r, target));
            assert(equal(velu(sum(q, generator, e), kernel.data(), kernel.size(), e), r));
            assert(equal(velu(sum(q, q, e), kernel.data(), kernel.size(), e), sum(r, r, target)));
            ++checked;
        }
    }
    std::cout << "p=" << p << " degree=" << ell << " source_a=" << e.a << " source_b=" << e.b
              << " order=" << order << " target_a=" << target.a << " target_b=" << target.b
              << " generator_x=" << generator.x << " generator_y=" << generator.y
              << " verified_points=" << checked << '\n';
}
