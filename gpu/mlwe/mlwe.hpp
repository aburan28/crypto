#pragma once
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace mlwe {
constexpr int N = 256;
constexpr int Q = 3329;

// Row-major layout for one instance: A[(row*k+col)*N+i],
// s[col*N+i], e[row*N+i], t[row*N+i]. Batches are contiguous.
inline std::size_t matrix_size(int k, int batch) { return std::size_t(batch)*k*k*N; }
inline std::size_t vector_size(int k, int batch) { return std::size_t(batch)*k*N; }
inline int canonical(std::int64_t x) {
    x %= Q;
    return static_cast<int>(x < 0 ? x + Q : x);
}

// Independent, simple reference: iterate over product terms and wrap x^256=-1.
inline std::vector<int> cpu_reference(const std::vector<int>& A,
                                      const std::vector<int>& s,
                                      const std::vector<int>& e,
                                      int k, int batch) {
    if (k < 2 || k > 4 || batch < 1 || A.size() != matrix_size(k,batch) ||
        s.size() != vector_size(k,batch) || e.size() != vector_size(k,batch))
        throw std::invalid_argument("invalid MLWE dimensions");
    std::vector<int> t(e.size());
    for (int b=0; b<batch; ++b) {
        for (int row=0; row<k; ++row) {
            std::int64_t accum[N]{};
            for (int col=0; col<k; ++col) {
                const int* a = A.data() + (std::size_t(b)*k*k + row*k + col)*N;
                const int* secret = s.data() + (std::size_t(b)*k + col)*N;
                for (int i=0; i<N; ++i)
                    for (int j=0; j<N; ++j) {
                        const std::int64_t term = std::int64_t(a[i])*secret[j];
                        const int degree = i+j;
                        accum[degree < N ? degree : degree-N] += degree < N ? term : -term;
                    }
            }
            for (int i=0; i<N; ++i) {
                const std::size_t off = (std::size_t(b)*k + row)*N+i;
                t[off] = canonical(accum[i]+e[off]);
            }
        }
    }
    return t;
}
} // namespace mlwe
