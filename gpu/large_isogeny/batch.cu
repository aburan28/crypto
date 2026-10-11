#include "velu.hpp"
#include <cuda_runtime.h>
#include <cstdint>

// Kernel points and inputs have already been certified on the host. Each
// thread evaluates the full degree-ell map at one independent input point.
__global__ void velu_batch(const iso::Point *inputs, iso::Point *outputs, uint64_t n,
                           const iso::Point *kernel, uint64_t kernel_count, iso::Curve source) {
    const uint64_t i = uint64_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < n) outputs[i] = iso::velu(inputs[i], kernel, kernel_count, source);
}

extern "C" cudaError_t iso_velu_batch(const iso::Point *device_inputs, iso::Point *device_outputs,
                                       uint64_t count, const iso::Point *device_kernel,
                                       uint64_t kernel_count, iso::Curve source, cudaStream_t stream) {
    if (!count) return cudaSuccess;
    if (!device_inputs || !device_outputs || !device_kernel || !kernel_count ||
        source.p < 5 || source.p >= (1ULL << 31)) return cudaErrorInvalidValue;
    velu_batch<<<(count + 127) / 128, 128, 0, stream>>>(device_inputs, device_outputs,
                                                       count, device_kernel, kernel_count, source);
    return cudaGetLastError();
}
