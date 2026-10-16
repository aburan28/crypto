#include "velu.hpp"
#include <cuda_runtime.h>
#include <cassert>
#include <iostream>
#include <vector>

extern "C" cudaError_t iso_velu_batch(const iso::Point *, iso::Point *, uint64_t,
                                       const iso::Point *, uint64_t, iso::Curve, cudaStream_t);

#define CUDA_OK(call) do { cudaError_t err = (call); if (err != cudaSuccess) { \
    std::cerr << #call << ": " << cudaGetErrorString(err) << '\n'; return 1; } } while (0)

int main() {
    const iso::Curve e{1009, 1, 7};
    const iso::Point g{573, 570, false};
    const iso::Curve target{1009, 753, 118};
    std::vector<iso::Point> kernel, inputs, outputs(1024);
    for (auto q = g; !q.infinity; q = iso::sum(q, g, e)) kernel.push_back(q);
    assert(kernel.size() == 100);
    for (uint64_t x = 0; x < e.p && inputs.size() < outputs.size(); ++x) {
        for (uint64_t y = 0; y < e.p && inputs.size() < outputs.size(); ++y) {
            iso::Point q{x, y, false};
            if (iso::on_curve(q, e)) inputs.push_back(q);
        }
    }
    assert(inputs.size() == outputs.size());
    iso::Point *di = nullptr, *do_ = nullptr, *dk = nullptr;
    CUDA_OK(cudaMalloc(&di, inputs.size() * sizeof(iso::Point)));
    CUDA_OK(cudaMalloc(&do_, outputs.size() * sizeof(iso::Point)));
    CUDA_OK(cudaMalloc(&dk, kernel.size() * sizeof(iso::Point)));
    CUDA_OK(cudaMemcpy(di, inputs.data(), inputs.size() * sizeof(iso::Point), cudaMemcpyHostToDevice));
    CUDA_OK(cudaMemcpy(dk, kernel.data(), kernel.size() * sizeof(iso::Point), cudaMemcpyHostToDevice));
    CUDA_OK(iso_velu_batch(di, do_, inputs.size(), dk, kernel.size(), e, nullptr));
    CUDA_OK(cudaDeviceSynchronize());
    CUDA_OK(cudaMemcpy(outputs.data(), do_, outputs.size() * sizeof(iso::Point), cudaMemcpyDeviceToHost));
    for (size_t i = 0; i < inputs.size(); ++i) {
        const auto expected = iso::velu(inputs[i], kernel.data(), kernel.size(), e);
        assert(iso::equal(outputs[i], expected));
        assert(iso::on_curve(outputs[i], target));
    }
    CUDA_OK(cudaFree(di)); CUDA_OK(cudaFree(do_)); CUDA_OK(cudaFree(dk));
    std::cout << "degree=101 gpu_cpu_equal=" << inputs.size() << '\n';
}
