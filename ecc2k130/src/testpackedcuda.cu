#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "../include/packed131.h"

using eccPacked131::P131;
static void checked(cudaError_t status) {
    if (status != cudaSuccess) {
        fprintf(stderr, "CUDA arithmetic test: %s\n", cudaGetErrorString(status));
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

struct RawPolynomial { uint32_t v[9]; };

static P131 reduceReference(RawPolynomial h) {
    // Long division, independent of the generated reciprocal reduction.
    const int terms[]={0,2,3,64,66,67,96,98,99,112,114,115,120,122,123,124,128,130,131};
    for (int degree=260;degree>=131;--degree) if ((h.v[degree/32]>>(degree%32))&1u) {
        for (int term:terms) {
            int bit=degree-131+term;
            h.v[bit/32]^=1u<<(bit%32);
        }
    }
    return P131{{h.v[0],h.v[1],h.v[2],h.v[3],h.v[4]&7u}};
}

static P131 multiplyReference(P131 a, P131 b) {
    RawPolynomial h{};
    for (int i=0;i<131;i++) if ((a.v[i/32]>>(i%32))&1u)
        for (int j=0;j<131;j++) if ((b.v[j/32]>>(j%32))&1u)
            h.v[(i+j)/32]^=1u<<((i+j)%32);
    return reduceReference(h);
}

__global__ __launch_bounds__(ECC_THREADS, ECC_MINBLOCKS)
void reductionProbe(const RawPolynomial *input, P131 *output, int n) {
    int i=blockIdx.x*blockDim.x+threadIdx.x;
    if (i<n) output[i]=eccPacked131::reducePolynomial131(input[i].v);
}

__global__ __launch_bounds__(ECC_THREADS, ECC_MINBLOCKS)
void multiplicationProbe(const P131 *a, const P131 *b, P131 *output, int n) {
    int i=blockIdx.x*blockDim.x+threadIdx.x;
    if (i<n) output[i]=eccPacked131::mulPolynomial131(a[i],b[i]);
}

static bool same(P131 a, P131 b) {
    for (int i=0;i<5;i++) if (a.v[i]!=b.v[i]) return false;
    return true;
}

static bool polynomialChecks() {
    uint32_t state=0x261131u;
    auto random=[&]() { state^=state<<13;state^=state>>17;state^=state<<5;return state; };
    std::vector<RawPolynomial> raw;
    for (int bit=0;bit<261;bit++) {
        RawPolynomial h{};h.v[bit/32]=1u<<(bit%32);raw.push_back(h);
    }
    for (int i=0;i<1000;i++) {
        RawPolynomial h;for (int j=0;j<9;j++) h.v[j]=random();h.v[8]&=31;raw.push_back(h);
    }
    int n=int(raw.size());std::vector<P131> output(n);
    RawPolynomial *deviceRaw;P131 *deviceOutput;
    checked(cudaMalloc(&deviceRaw,n*sizeof(RawPolynomial)));
    checked(cudaMalloc(&deviceOutput,n*sizeof(P131)));
    checked(cudaMemcpy(deviceRaw,raw.data(),n*sizeof(RawPolynomial),cudaMemcpyHostToDevice));
    reductionProbe<<<(n+ECC_THREADS-1)/ECC_THREADS,ECC_THREADS>>>(deviceRaw,deviceOutput,n);
    checked(cudaGetLastError());checked(cudaDeviceSynchronize());
    checked(cudaMemcpy(output.data(),deviceOutput,n*sizeof(P131),cudaMemcpyDeviceToHost));
    checked(cudaFree(deviceRaw));checked(cudaFree(deviceOutput));
    for (int i=0;i<n;i++) if (!same(output[i],reduceReference(raw[i]))) {
        fprintf(stderr,"GPU polynomial reduction mismatch at %d\n",i);return false;
    }
    printf("PASS: %d GPU polynomial reductions against long division\n",n);

    std::vector<P131> a,b;
    for (int i=0;i<131;i++) for (int j=0;j<131;j++) {
        P131 x{},y{};x.v[i/32]=1u<<(i%32);y.v[j/32]=1u<<(j%32);a.push_back(x);b.push_back(y);
    }
    for (int i=0;i<1024;i++) {
        P131 x,y;for (int j=0;j<5;j++) {x.v[j]=random();y.v[j]=random();}
        x.v[4]&=7;y.v[4]&=7;a.push_back(x);b.push_back(y);
    }
    n=int(a.size());output.resize(n);P131 *deviceA,*deviceB;
    checked(cudaMalloc(&deviceA,n*sizeof(P131)));checked(cudaMalloc(&deviceB,n*sizeof(P131)));
    checked(cudaMalloc(&deviceOutput,n*sizeof(P131)));
    checked(cudaMemcpy(deviceA,a.data(),n*sizeof(P131),cudaMemcpyHostToDevice));
    checked(cudaMemcpy(deviceB,b.data(),n*sizeof(P131),cudaMemcpyHostToDevice));
    multiplicationProbe<<<(n+ECC_THREADS-1)/ECC_THREADS,ECC_THREADS>>>(deviceA,deviceB,deviceOutput,n);
    checked(cudaGetLastError());checked(cudaDeviceSynchronize());
    checked(cudaMemcpy(output.data(),deviceOutput,n*sizeof(P131),cudaMemcpyDeviceToHost));
    checked(cudaFree(deviceA));checked(cudaFree(deviceB));checked(cudaFree(deviceOutput));
    for (int i=0;i<n;i++) if (!same(output[i],multiplyReference(a[i],b[i]))) {
        fprintf(stderr,"GPU polynomial multiplication mismatch at %d\n",i);return false;
    }
    printf("PASS: %d GPU polynomial products, including all 17161 basis pairs\n",n);
    return true;
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
    return polynomialChecks()?0:1;
}
