#include "mlwe.hpp"
#include <cuda_runtime.h>
#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>

#define CUDA_CHECK(call) do { cudaError_t err=(call); if(err!=cudaSuccess) { \
    std::cerr << #call << ": " << cudaGetErrorString(err) << "\n"; std::exit(2); \
} } while(0)

// One thread computes one coefficient. Correct reference implementation, with
// 64-bit sums to support any canonical residue in s, as well as short secrets.
__global__ void module_lwe_kernel(const int* __restrict__ A,
                                  const int* __restrict__ s,
                                  const int* __restrict__ e,
                                  int* __restrict__ t,
                                  int k) {
    const int slot = blockIdx.x*blockDim.x + threadIdx.x;
    if (slot >= k*mlwe::N) return;
    const int b=blockIdx.y;
    const int row=slot/mlwe::N;
    const int degree=slot%mlwe::N;
    std::int64_t sum=e[(std::size_t(b)*k+row)*mlwe::N+degree];
    for (int col=0;col<k;++col) {
        const int* a=A+(std::size_t(b)*k*k+row*k+col)*mlwe::N;
        const int* secret=s+(std::size_t(b)*k+col)*mlwe::N;
        for (int i=0;i<mlwe::N;++i) {
            const int j=degree>=i ? degree-i : degree-i+mlwe::N;
            const std::int64_t term=std::int64_t(a[i])*secret[j];
            sum += degree>=i ? term : -term;
        }
    }
    int v=static_cast<int>(sum%mlwe::Q);
    t[(std::size_t(b)*k+row)*mlwe::N+degree]=v<0 ? v+mlwe::Q : v;
}

void gpu_product(const std::vector<int>& A,const std::vector<int>& s,
                 const std::vector<int>& e,std::vector<int>& t,int k,int batch) {
    if (k<2 || k>4 || batch<1 || batch>65535 ||
        A.size()!=mlwe::matrix_size(k,batch) || s.size()!=mlwe::vector_size(k,batch)
        || e.size()!=s.size() || t.size()!=s.size())
        throw std::invalid_argument("invalid CUDA MLWE dimensions");
    int *da=nullptr,*ds=nullptr,*de=nullptr,*dt=nullptr;
    CUDA_CHECK(cudaMalloc(&da,A.size()*sizeof(int)));
    CUDA_CHECK(cudaMalloc(&ds,s.size()*sizeof(int)));
    CUDA_CHECK(cudaMalloc(&de,e.size()*sizeof(int)));
    CUDA_CHECK(cudaMalloc(&dt,t.size()*sizeof(int)));
    CUDA_CHECK(cudaMemcpy(da,A.data(),A.size()*sizeof(int),cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(ds,s.data(),s.size()*sizeof(int),cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(de,e.data(),e.size()*sizeof(int),cudaMemcpyHostToDevice));
    constexpr int threads=128;
    const dim3 blocks((k*mlwe::N+threads-1)/threads,batch);
    module_lwe_kernel<<<blocks,threads>>>(da,ds,de,dt,k);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaMemcpy(t.data(),dt,t.size()*sizeof(int),cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaFree(da)); CUDA_CHECK(cudaFree(ds));
    CUDA_CHECK(cudaFree(de)); CUDA_CHECK(cudaFree(dt));
}

int main() {
    int count=0;
    CUDA_CHECK(cudaGetDeviceCount(&count));
    if (!count) { std::cerr << "No CUDA GPU available\n"; return 2; }
    for (int k : {2,3,4}) {
        constexpr int batch=8;
        std::vector<int> A(mlwe::matrix_size(k,batch));
        std::vector<int> s(mlwe::vector_size(k,batch));
        std::vector<int> e(s.size()),t(s.size());
        std::mt19937 rng(314159+k); // Deterministic fixtures only.
        std::uniform_int_distribution<int> u(0,mlwe::Q-1),small(-3,3);
        for (int& x:A) x=u(rng);
        for (int& x:s) x=small(rng);
        for (int& x:e) x=small(rng);
        auto expected=mlwe::cpu_reference(A,s,e,k,batch);
        gpu_product(A,s,e,t,k,batch);
        if (t!=expected) {
            for (std::size_t i=0;i<t.size();++i)
                if (t[i]!=expected[i]) {
                    std::cerr << "FAIL k=" << k << " index=" << i
                              << " gpu=" << t[i] << " cpu=" << expected[i] << "\n";
                    break;
                }
            return 1;
        }
        std::cout << "PASS CUDA vs CPU k=" << k << " batch=" << batch << "\n";
    }
}
