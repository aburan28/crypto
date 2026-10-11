#include "spmv_io.hpp"

#include <cuda_runtime.h>
#include <iostream>

static void check_cuda(cudaError_t status, const char *operation) {
    if (status != cudaSuccess)
        throw std::runtime_error(
            std::string(operation) + ": " + cudaGetErrorString(status));
}

__device__ static std::uint64_t add_mod(
    std::uint64_t a, std::uint64_t b, std::uint64_t modulus) {
    const std::uint64_t sum = a + b; // a,b < 2^63, so this cannot overflow.
    return sum >= modulus ? sum - modulus : sum;
}

__device__ static std::uint64_t multiply_mod(
    std::uint64_t a, std::uint64_t b, std::uint64_t modulus) {
    std::uint64_t result = 0;
    while (b != 0) {
        if (b & 1) result = add_mod(result, a, modulus);
        b >>= 1;
        if (b != 0) a = add_mod(a, a, modulus);
    }
    return result;
}

__global__ static void spmv_kernel(
    const std::uint64_t *row_ptr,
    const std::uint64_t *column_index,
    const std::uint64_t *coefficient,
    const std::uint64_t *x,
    std::uint64_t *y,
    std::uint64_t rows,
    std::uint64_t lanes,
    std::uint64_t modulus) {
    const std::uint64_t task = blockIdx.x * blockDim.x + threadIdx.x;
    if (task >= rows * lanes) return;
    const std::uint64_t row = task / lanes;
    const std::uint64_t lane = task % lanes;
    std::uint64_t sum = 0;
    for (std::uint64_t at = row_ptr[row]; at < row_ptr[row + 1]; ++at) {
        const auto column = column_index[at];
        sum = add_mod(
            sum,
            multiply_mod(coefficient[at] % modulus, x[column * lanes + lane], modulus),
            modulus);
    }
    y[task] = sum;
}

template <class T>
static T *copy_to_device(const std::vector<T> &host) {
    T *device = nullptr;
    check_cuda(cudaMalloc(&device, host.size() * sizeof(T)), "cudaMalloc");
    check_cuda(
        cudaMemcpy(device, host.data(), host.size() * sizeof(T), cudaMemcpyHostToDevice),
        "cudaMemcpy host-to-device");
    return device;
}

int main(int argc, char **argv) {
    try {
        const auto p = read_spmv1(input_path(argc, argv));
        auto *row_ptr = copy_to_device(p.row_ptr);
        auto *column_index = copy_to_device(p.column_index);
        auto *coefficient = copy_to_device(p.coefficient);
        auto *x = copy_to_device(p.x);
        std::uint64_t *y = nullptr;
        check_cuda(
            cudaMalloc(&y, p.rows * p.lanes * sizeof(std::uint64_t)),
            "cudaMalloc output");
        const std::uint64_t tasks = p.rows * p.lanes;
        spmv_kernel<<<(tasks + 255) / 256, 256>>>(
            row_ptr, column_index, coefficient, x, y, p.rows, p.lanes, p.modulus);
        check_cuda(cudaGetLastError(), "launch SPMV1 kernel");
        check_cuda(cudaDeviceSynchronize(), "run SPMV1 kernel");
        std::vector<std::uint64_t> output(tasks);
        check_cuda(
            cudaMemcpy(output.data(), y, tasks * sizeof(std::uint64_t), cudaMemcpyDeviceToHost),
            "cudaMemcpy device-to-host");
        cudaFree(row_ptr);
        cudaFree(column_index);
        cudaFree(coefficient);
        cudaFree(x);
        cudaFree(y);
        print_spmv1(p, output);
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
