// The table walk's cycle rule (tablewalk.h, eccTagFruitless) on its own, with
// no curve: which runs of steps it lets a walk complete.
//
// A run of step tags returns a walk to a point it has left exactly when, per
// branch h, the signed sum of sigma^k over the branch's steps is O -- in pairs
// (a tag and its negation), or through sigma^2 + sigma + 2 = 0 and
// sigma^3 + sigma - 2 = 0 (WALK-CONSTANT.md section 5).  This walks each such
// run through the rule exactly as the walk would, history and all, and checks
// that the rule refuses a step of every run it must break:
//
//   - every tau-relation 4-cycle (24 orders), at every phase including the
//     wrap mod m;
//   - every 6-step run that returns through Frobenius without a pair (the
//     nine families of benchmarks/walk-constant, random orders and phases);
//   - every pairwise cycle of 2 to 10 steps (all perfect matchings);
//
// and lets through the shortest pairwise cycle it is not meant to see, the
// 12-step a..f -a..-f, so the boundary is where the note says.  Then: an
// empty history never fires, random tags rarely fire without a pair, and the
// rule's answer is unchanged when every tag is conjugated by sigma^c or
// negated -- the rule has to be a class function, as the step is.
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>
#include "../include/tablewalk.h"

// eccTagFruitless with the cycle tables for degree m, as TableWalkConsts
// builds them.
static bool fruitless(unsigned t, unsigned long long hist, int m) {
    static std::vector<uint16_t> rpow;
    static int built = 0;
    if (built != m) {
        rpow.assign(m + 1, 0);
        for (int k = 0; k < m; ++k) rpow[k] = uint16_t(eccCyclePow(k));
        built = m;
    }
    const EccCycleWindow w = eccCycleWindow(eccTagK(t), hist, rpow.data(), m, eccCyclePow(m));
    return eccTagFruitless(t, w, rpow.data(), m, eccCyclePow(m));
}

static int failures = 0;
static void check(bool ok, const char *what) {
    if (!ok) {
        ++failures;
        std::printf("FAIL: %s\n", what);
    }
}

// Walk a run of tags through the rule, starting from `start` history.  Returns
// the index of the first step the rule refuses, or -1 if it takes every step.
static int firstRefusal(const std::vector<unsigned> &tags, unsigned long long start, int m) {
    unsigned long long hist = start;
    for (size_t i = 0; i < tags.size(); ++i) {
        if (fruitless(tags[i], hist, m)) return int(i);
        hist = eccHistPush(hist, tags[i]);
    }
    return -1;
}

static unsigned conj(unsigned t, int c, int m, bool negate) {
    if (t == ECC_TAG_NONE) return t;
    const int k = (eccTagK(t) + c) % m;
    return eccTag(eccTagH(t), k, eccTagEps(t) ^ (negate ? 1 : 0));
}

static void matchings(std::vector<int> &pos, std::vector<std::pair<int, int>> &cur,
                      std::vector<std::vector<std::pair<int, int>>> &out) {
    if (pos.empty()) { out.push_back(cur); return; }
    const int a = pos[0];
    for (size_t i = 1; i < pos.size(); ++i) {
        std::vector<int> rest;
        for (size_t j = 1; j < pos.size(); ++j) if (j != i) rest.push_back(pos[j]);
        cur.push_back({a, pos[i]});
        matchings(rest, cur, out);
        cur.pop_back();
    }
}

int main() {
    std::mt19937_64 rng(20260924);
    const int H = 8;
    for (int m : {23, 41, 131}) {
        std::uniform_int_distribution<int> hDist(0, H - 1), kDist(0, m - 1), eDist(0, 1);
        auto randomTag = [&]() { return eccTag(hDist(rng), kDist(rng), eDist(rng)); };
        // A random history of real tags, as a walk mid-trail has.
        auto randomHistory = [&]() {
            unsigned long long h = ECC_HIST_EMPTY;
            for (int i = 0; i < ECC_HIST_DEPTH; ++i) h = eccHistPush(h, randomTag());
            return h;
        };

        // 1. Every tau-relation 4-cycle, every order, phases across the wrap.
        int tauRuns = 0, tauBroken = 0;
        for (int rel = 0; rel < 2; ++rel) {
            for (int trial = 0; trial < 400; ++trial) {
                const int h = hDist(rng), e = eDist(rng);
                const int k = trial < 8 ? m - 1 - trial % 4 : kDist(rng);   // wrap cases first
                const unsigned r = eccTag(h, k, e);
                std::vector<unsigned> base = {r, r};
                if (rel == 0) {
                    base.push_back(eccTagAdvanceK(r, 1, m));
                    base.push_back(eccTagAdvanceK(r, 2, m));
                } else {
                    base.push_back(eccTagAdvanceK(r, 1, m) ^ ECC_TAG_EPS);
                    base.push_back(eccTagAdvanceK(r, 3, m) ^ ECC_TAG_EPS);
                }
                std::vector<int> idx = {0, 1, 2, 3};
                do {
                    std::vector<unsigned> run;
                    for (int i : idx) run.push_back(base[i]);
                    // Two laps: a walk in the cycle would repeat it.
                    std::vector<unsigned> laps = run;
                    laps.insert(laps.end(), run.begin(), run.end());
                    ++tauRuns;
                    const int at = firstRefusal(laps, randomHistory(), m);
                    tauBroken += at >= 0 && at < 4;
                } while (std::next_permutation(idx.begin(), idx.end()));
            }
        }
        std::printf("m = %d: tau-relation 4-cycles refused within the first lap: %d / %d\n", m, tauBroken,
                    tauRuns);
        check(tauBroken == tauRuns, "every tau-relation 4-cycle is refused within its first lap");

        // 1b. The 6-step runs that return through Frobenius with no pair:
        // (sign, offset) terms of benchmarks/walk-constant/fruitless_patterns.
        static const int six[9][6][2] = {
            {{1, 0}, {1, 0}, {1, 1}, {-1, 2}, {1, 3}, {1, 5}},
            {{1, 0}, {1, 0}, {-1, 1}, {-1, 2}, {-1, 2}, {1, 5}},
            {{1, 0}, {1, 0}, {1, 0}, {1, 0}, {-1, 2}, {1, 5}},
            {{1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 2}, {-1, 3}},
            {{1, 0}, {1, 0}, {1, 1}, {-1, 2}, {-1, 3}, {-1, 4}},
            {{1, 0}, {1, 0}, {-1, 1}, {1, 3}, {-1, 4}, {-1, 6}},
            {{1, 0}, {1, 0}, {-1, 1}, {1, 2}, {1, 2}, {1, 4}},
            {{1, 0}, {1, 0}, {-1, 1}, {1, 3}, {1, 4}, {1, 5}},
            {{1, 0}, {1, 0}, {1, 1}, {1, 1}, {1, 1}, {-1, 4}},
        };
        int sixRuns = 0, sixBroken = 0;
        for (int f = 0; f < 9; ++f)
            for (int trial = 0; trial < 200; ++trial) {
                const int h = hDist(rng), e = eDist(rng);
                const int k = trial < 8 ? m - 1 - trial % 6 : kDist(rng);
                std::vector<unsigned> run;
                for (int j = 0; j < 6; ++j)
                    run.push_back(eccTag(h, (k + six[f][j][1]) % m, e ^ (six[f][j][0] < 0 ? 1 : 0)));
                std::shuffle(run.begin(), run.end(), rng);
                std::vector<unsigned> laps = run;
                laps.insert(laps.end(), run.begin(), run.end());
                ++sixRuns;
                const int at = firstRefusal(laps, randomHistory(), m);
                sixBroken += at >= 0 && at < 6;
            }
        std::printf("m = %d: 6-step Frobenius runs refused within the first lap: %d / %d\n", m, sixBroken,
                    sixRuns);
        check(sixBroken == sixRuns, "every 6-step Frobenius run is refused within its first lap");

        // 2. Pairwise cycles of 2..10 steps are broken; the 12-step one is not.
        for (int L = 2; L <= 12; L += 2) {
            std::vector<int> pos;
            for (int i = 0; i < L; ++i) pos.push_back(i);
            std::vector<std::vector<std::pair<int, int>>> all;
            std::vector<std::pair<int, int>> cur;
            matchings(pos, cur, all);
            int broken = 0, total = 0;
            for (const auto &mt : all) {
                if (L == 12) {
                    // Only the all-distance-6 matching is the boundary case.
                    bool far = true;
                    for (auto pr : mt) far &= pr.second - pr.first == 6;
                    if (!far) continue;
                }
                if (L == 10 && total >= 2000) break;
                for (int trial = 0; trial < (L <= 10 ? 20 : 200); ++trial) {
                    std::vector<unsigned> run(L);
                    for (auto pr : mt) {
                        // Distinct branches per pair, so no tau-relation and no
                        // accidental extra cancellation.
                        const unsigned a = eccTag((pr.first * 3 + trial) % H, kDist(rng), eDist(rng));
                        run[pr.first] = a;
                        run[pr.second] = a ^ ECC_TAG_EPS;
                    }
                    std::vector<unsigned> laps = run;
                    laps.insert(laps.end(), run.begin(), run.end());
                    ++total;
                    const int at = firstRefusal(laps, randomHistory(), m);
                    broken += at >= 0 && at < L;
                }
            }
            if (L <= 10) {
                std::printf("m = %d: pairwise %d-cycles refused within the first lap: %d / %d\n", m, L, broken,
                            total);
                check(broken == total, "every pairwise cycle of at most 10 steps is refused");
            } else {
                std::printf("m = %d: the 12-step a..f,-a..-f cycle refused in %d / %d (expected rarely: "
                            "only by an accidental match with the history)\n", m, broken, total);
                check(broken < total / 10, "the 12-step pairwise cycle is outside the rule, as documented");
            }
        }

        // 3. An empty history never fires; random tags fire only on a pair.
        bool emptyOk = true;
        for (int h = 0; h < H; ++h)
            for (int k = 0; k < m; ++k)
                for (int e = 0; e < 2; ++e) emptyOk &= !fruitless(eccTag(h, k, e), ECC_HIST_EMPTY, m);
        check(emptyOk, "an empty history never fires");
        int spurious = 0;
        for (int i = 0; i < 200000; ++i) {
            const unsigned t = randomTag();
            const unsigned long long hist = randomHistory();
            bool pair = false;
            for (int j = 0; j < ECC_HIST_DEPTH; ++j)
                pair |= eccTagNegates(t, unsigned(hist >> (ECC_HIST_SLOT * j)) & ECC_TAG_MASK);
            // Without a pair, only a Frobenius relation (about 8^-3 (2m)^-3
            // per window) or the 2^-16 of the Z/2^16 image can fire.
            spurious += !pair && fruitless(t, hist, m);
        }
        std::printf("m = %d: 200,000 random steps fired without a pair %d times\n", m, spurious);
        check(spurious <= 40, "random tags fire without a pair only at the Z/2^16 rate");

        // 4. Class function: conjugating every tag by sigma^c and/or negating
        // it leaves the answer unchanged, on random and on firing histories.
        int covBad = 0, covFired = 0;
        for (int i = 0; i < 100000; ++i) {
            unsigned hs[ECC_HIST_DEPTH];
            for (int j = 0; j < ECC_HIST_DEPTH; ++j) hs[j] = randomTag();
            unsigned t = randomTag();
            if (i % 3 == 1) t = hs[1 + i % 3] ^ ECC_TAG_EPS;                   // pairwise fire
            if (i % 3 == 2) {                                                  // tau fire
                hs[0] = t;
                hs[1] = eccTagAdvanceK(t, 1, m);
                hs[2] = eccTagAdvanceK(t, 2, m);
            }
            unsigned long long hist = ECC_HIST_EMPTY;
            for (int j = ECC_HIST_DEPTH - 1; j >= 0; --j) hist = eccHistPush(hist, hs[j]);
            const bool f = fruitless(t, hist, m);
            covFired += f;
            const int c = 1 + int(rng() % (m - 1));
            for (int neg = 0; neg < 2; ++neg) {
                unsigned long long h2 = ECC_HIST_EMPTY;
                for (int j = ECC_HIST_DEPTH - 1; j >= 0; --j) h2 = eccHistPush(h2, conj(hs[j], c, m, neg));
                covBad += fruitless(conj(t, c, m, neg), h2, m) != f;
            }
        }
        std::printf("m = %d: covariance under sigma^c and negation: %d mismatches in 200,000 (%d fired)\n", m,
                    covBad, covFired);
        check(covBad == 0, "the rule is a class function");
        check(covFired >= 60000, "the covariance probe exercises firing histories");
    }
    std::printf("%s\n", failures ? "FAIL" : "PASS");
    return failures ? 1 : 0;
}
