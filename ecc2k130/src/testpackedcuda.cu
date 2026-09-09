#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "../include/packed131.h"

using eccPacked131::P131;
static void checked(cudaError_t status) {
    if (status != cudaSuccess) {
        fprintf(stderr, "CUDA Frobenius test: %s\n", cudaGetErrorString(status));
        std::exit(1);
    }
}

__global__ __launch_bounds__(ECC_THREADS, ECC_MINBLOCKS)
void frobeniusProbe(const P131 *input, const int *powers, P131 *output, int n) {
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    if (i < n) output[i] = eccPacked131::sigma131(input[i], powers[i]);
}

static P131 expected(P131 a, int power) {
    unsigned factor = 1;
    for (int i=0;i<power;i++) factor = (2*factor)%263;
    P131 out{};
    for (int i=0;i<131;i++) if ((a.v[i/32]>>(i%32))&1u) {
        unsigned j = ((i+1)*factor)%263;
        if (j>131) j=263-j;
        --j;
        out.v[j/32] ^= 1u<<(j%32);
    }
    return out;
}

int main() {
    const int selected[] = {0,1,2,3,4,5,6,7,8,9,10,16,32,65,130,131};
    std::vector<P131> input;
    std::vector<int> powers;
    for (int k:selected) for (int bit=0;bit<131;bit++) {
        P131 a{}; a.v[bit/32]=1u<<(bit%32);
        input.push_back(a); powers.push_back(k);
    }
    uint32_t state=0x131263u;
    for (int i=0;i<1024;i++) {
        P131 a;
        for (int j=0;j<5;j++) {
            state^=state<<13; state^=state>>17; state^=state<<5;
            a.v[j]=state;
        }
        a.v[4]&=7;
        input.push_back(a); powers.push_back(selected[i%16]);
    }
    int n=int(input.size());
    std::vector<P131> output(n);
    P131 *deviceInput, *deviceOutput; int *devicePowers;
    checked(cudaMalloc(&deviceInput, n*sizeof(P131)));
    checked(cudaMalloc(&deviceOutput, n*sizeof(P131)));
    checked(cudaMalloc(&devicePowers, n*sizeof(int)));
    checked(cudaMemcpy(deviceInput,input.data(),n*sizeof(P131),cudaMemcpyHostToDevice));
    checked(cudaMemcpy(devicePowers,powers.data(),n*sizeof(int),cudaMemcpyHostToDevice));
    frobeniusProbe<<<(n+ECC_THREADS-1)/ECC_THREADS,ECC_THREADS>>>(deviceInput,devicePowers,deviceOutput,n);
    checked(cudaGetLastError()); checked(cudaDeviceSynchronize());
    checked(cudaMemcpy(output.data(),deviceOutput,n*sizeof(P131),cudaMemcpyDeviceToHost));
    checked(cudaFree(deviceInput)); checked(cudaFree(deviceOutput)); checked(cudaFree(devicePowers));
    for (int i=0;i<n;i++) {
        P131 want=expected(input[i],powers[i]);
        for (int j=0;j<5;j++) if (output[i].v[j]!=want.v[j]) {
            fprintf(stderr,"GPU Frobenius mismatch: case %d, power %d, word %d\n",i,powers[i],j);
            return 1;
        }
    }
    printf("PASS: %d GPU Frobenius vectors, every field basis vector for all selected powers plus dense cases\n",n);
}
