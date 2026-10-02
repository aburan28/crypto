#include <cuda_runtime.h>

#include "../include/curveparams.h"
#include "../include/packed131.h"
#include "../include/packeddirectsigma131.cuh"

#include <cstdio>
#include <vector>

using P = eccPacked131::P131;

#define CUDA_CHECK(call) do { \
    const cudaError_t error_ = (call); \
    if (error_ != cudaSuccess) { \
        std::fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(error_)); \
        return 2; \
    } \
} while (0)

__global__ void directSigmaProbe(const P *input, const unsigned char *indices, P *output, int count) {
    extern __shared__ uint32_t table[];
    eccPacked131::initDirectSigmaShared131(table);
    for (int i = int(threadIdx.x); i < count; i += int(blockDim.x))
        output[i] = eccPacked131::directSigmaShared131(input[i], indices[i], table);
}

static bool same(const P &a, const P &b) {
    for (int word = 0; word < 5; ++word)
        if (a.v[word] != b.v[word]) return false;
    return (a.v[4] & ~7u) == 0;
}

static P oracle(P input, int index) {
    const P normal = eccPacked131::fromPolynomial131(input);
    return eccPacked131::toPolynomial131(
        eccPacked131::add131(normal, eccPacked131::sigma131(normal, index + 3)));
}

static uint32_t nextWord(uint32_t &state) {
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return state;
}

int main() {
    constexpr int densePerMap = 512;
    std::vector<P> input;
    std::vector<unsigned char> indices;
    std::vector<P> expected;
    uint32_t state = 0xd1ec7513u;
    for (int index = 0; index < 8; ++index) {
        for (int bit = 0; bit < 131; ++bit) {
            P value{{0, 0, 0, 0, 0}};
            value.v[bit / 32] = 1u << (bit % 32);
            input.push_back(value);
            indices.push_back(static_cast<unsigned char>(index));
            expected.push_back(oracle(value, index));
        }
        for (int test = 0; test < densePerMap; ++test) {
            P value{{nextWord(state), nextWord(state), nextWord(state), nextWord(state), nextWord(state) & 7u}};
            if (test == 0) value = P{{0, 0, 0, 0, 0}};
            if (test == 1) value = P{{~0u, ~0u, ~0u, ~0u, 7u}};
            input.push_back(value);
            indices.push_back(static_cast<unsigned char>(index));
            expected.push_back(oracle(value, index));
        }
    }

    P *deviceInput = nullptr, *deviceOutput = nullptr;
    unsigned char *deviceIndices = nullptr;
    const size_t fieldBytes = input.size() * sizeof(P);
    CUDA_CHECK(cudaMalloc(&deviceInput, fieldBytes));
    CUDA_CHECK(cudaMalloc(&deviceOutput, fieldBytes));
    CUDA_CHECK(cudaMalloc(&deviceIndices, indices.size()));
    CUDA_CHECK(cudaMemcpy(deviceInput, input.data(), fieldBytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(deviceIndices, indices.data(), indices.size(), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaFuncSetAttribute(directSigmaProbe, cudaFuncAttributeMaxDynamicSharedMemorySize,
                                    eccPacked131::DIRECT_SIGMA_SHARED_BYTES));
    directSigmaProbe<<<1, 512, eccPacked131::DIRECT_SIGMA_SHARED_BYTES>>>(
        deviceInput, deviceIndices, deviceOutput, int(input.size()));
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    std::vector<P> output(input.size());
    CUDA_CHECK(cudaMemcpy(output.data(), deviceOutput, fieldBytes, cudaMemcpyDeviceToHost));
    for (size_t i = 0; i < output.size(); ++i) {
        if (!same(output[i], expected[i])) {
            std::fprintf(stderr, "direct sigma mismatch at case %zu, j=%d\n", i, int(indices[i]) + 3);
            return 1;
        }
    }
    cudaFree(deviceInput);
    cudaFree(deviceOutput);
    cudaFree(deviceIndices);
    std::printf("PASS: %zu GPU direct-sigma cases (%d basis and %d dense)\n",
                output.size(), 8 * 131, 8 * densePerMap);
    return 0;
}
