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
//   - every pairwise cycle of 2 to 8 steps (all perfect matchings);
//
// and lets through the shortest pairwise cycle it is not meant to see, the
// 10-step a b c d e -a -b -c -d -e, so the boundary is where the note says.
// Then: an empty history never fires, random tags never trip the tau test,
// and the rule's answer is unchanged when every tag is conjugated by sigma^c
// or negated -- the rule has to be a class function, as the step is.
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>
#include "../include/tablewalk.h"

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
        if (eccTagFruitless(tags[i], hist, m)) return int(i);
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
            for (int i = 0; i < 4; ++i) h = eccHistPush(h, randomTag());
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

        // 2. Pairwise cycles of 2..8 steps are broken; the 10-step one is not.
        for (int L = 2; L <= 10; L += 2) {
            std::vector<int> pos;
            for (int i = 0; i < L; ++i) pos.push_back(i);
            std::vector<std::vector<std::pair<int, int>>> all;
            std::vector<std::pair<int, int>> cur;
            matchings(pos, cur, all);
            int broken = 0, total = 0;
            for (const auto &mt : all) {
                if (L == 10) {
                    // Only the all-distance-5 matching is the boundary case.
                    bool far = true;
                    for (auto pr : mt) far &= pr.second - pr.first == 5;
                    if (!far) continue;
                }
                for (int trial = 0; trial < (L <= 8 ? 20 : 200); ++trial) {
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
            if (L <= 8) {
                std::printf("m = %d: pairwise %d-cycles refused within the first lap: %d / %d\n", m, L, broken,
                            total);
                check(broken == total, "every pairwise cycle of at most 8 steps is refused");
            } else {
                std::printf("m = %d: the 10-step a..e,-a..-e cycle refused in %d / %d (expected rarely: "
                            "only by an accidental match with the history)\n", m, broken, total);
                check(broken < total / 10, "the 10-step pairwise cycle is outside the rule, as documented");
            }
        }

        // 3. An empty history never fires; random tags never trip the tau test.
        bool emptyOk = true;
        for (int h = 0; h < H; ++h)
            for (int k = 0; k < m; ++k)
                for (int e = 0; e < 2; ++e) emptyOk &= !eccTagFruitless(eccTag(h, k, e), ECC_HIST_EMPTY, m);
        check(emptyOk, "an empty history never fires");
        int tauRandom = 0;
        for (int i = 0; i < 200000; ++i) {
            const unsigned t = randomTag(), a = randomTag(), b = randomTag(), c = randomTag();
            // Count only the tau test on tags that are not a relation by
            // construction: a random quadruple is one with probability
            // about 24 * 8^-3 * (2m)^-3, far below one in 200,000.
            tauRandom += eccTauRelation(t, a, b, c, m);
        }
        std::printf("m = %d: tau test on 200,000 random quadruples fired %d times\n", m, tauRandom);
        check(tauRandom <= 1, "random tags do not trip the tau test");

        // 4. Class function: conjugating every tag by sigma^c and/or negating
        // it leaves the answer unchanged, on random and on firing histories.
        int covBad = 0, covFired = 0;
        for (int i = 0; i < 100000; ++i) {
            unsigned hs[4];
            for (int j = 0; j < 4; ++j) hs[j] = randomTag();
            unsigned t = randomTag();
            if (i % 3 == 1) t = hs[1 + i % 3] ^ ECC_TAG_EPS;                   // pairwise fire
            if (i % 3 == 2) {                                                  // tau fire
                hs[0] = t;
                hs[1] = eccTagAdvanceK(t, 1, m);
                hs[2] = eccTagAdvanceK(t, 2, m);
            }
            unsigned long long hist = ECC_HIST_EMPTY;
            for (int j = 3; j >= 0; --j) hist = eccHistPush(hist, hs[j]);
            const bool f = eccTagFruitless(t, hist, m);
            covFired += f;
            const int c = 1 + int(rng() % (m - 1));
            for (int neg = 0; neg < 2; ++neg) {
                unsigned long long h2 = ECC_HIST_EMPTY;
                for (int j = 3; j >= 0; --j) h2 = eccHistPush(h2, conj(hs[j], c, m, neg));
                covBad += eccTagFruitless(conj(t, c, m, neg), h2, m) != f;
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
