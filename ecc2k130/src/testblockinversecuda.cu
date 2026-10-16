#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include "../include/packedblockinverse131.cuh"
using eccPacked131::P131;
static void checked(cudaError_t status) {
    if (status!=cudaSuccess) { fprintf(stderr,"CUDA block inverse test: %s\n",cudaGetErrorString(status));std::exit(1); }
}
template<bool PairZero>
__global__ __launch_bounds__(256,2)
void probe(const P131 *input,P131 *output,int n,int pitch,int guard) {
    __shared__ uint32_t tree[eccPacked131::blockInverseWords131 + (PairZero ? eccPacked131::blockInverseFlagWords131 : 0)];
    const int tid=blockIdx.x*blockDim.x+threadIdx.x;
    P131 a=tid<n?input[tid]:P131{{1,0,0,0,0}};
    P131 inverse=eccPacked131::blockInverse131<PairZero>(a,tree);
    if (tid<n) output[guard+tid]=inverse;
    P131 roundtrip=eccPacked131::blockInverse131<PairZero>(inverse,tree);
    if (tid<n) output[pitch+guard+tid]=roundtrip;
}
static P131 multiplyReference(P131 a,P131 b) {
    // Independent bit product and long division over the field polynomial.
    uint32_t h[9]={};
    for (int i=0;i<131;i++) if ((a.v[i/32]>>(i%32))&1u)
        for (int j=0;j<131;j++) if ((b.v[j/32]>>(j%32))&1u)
            h[(i+j)/32]^=1u<<((i+j)%32);
    const int terms[]={0,2,3,64,66,67,96,98,99,112,114,115,120,122,123,124,128,130,131};
    for (int degree=260;degree>=131;--degree) if ((h[degree/32]>>(degree%32))&1u)
        for (int term:terms) { int bit=degree-131+term;h[bit/32]^=1u<<(bit%32); }
    return P131{{h[0],h[1],h[2],h[3],h[4]}};
}
static bool same(P131 a,P131 b) { return std::memcmp(&a,&b,sizeof(a))==0; }
int main() {
    const int sizes[]={0,1,2,31,32,33,63,64,65,127,128,129,255,256,257,511,512,513,1025};
    const P131 one{{1,0,0,0,0}},zero{};
    P131 canary;std::memset(&canary,0xa5,sizeof(canary));
    uint32_t state=0x131263a5u;
    auto random=[&]() { state^=state<<13;state^=state>>17;state^=state<<5;return state; };
    size_t values=0,scenarios=0,zeros=0;
    for (int paired=0;paired<2;++paired) for (int pattern=0;pattern<4;++pattern) for (int n:sizes) {
        const int guard=31,pitch=((n+255)/256+1)*256+2*guard;
        std::vector<P131> input(n?n:1),output(2*pitch,canary);
        for (int i=0;i<n;++i) {
            P131 a;
            for (int w=0;w<5;++w) a.v[w]=random();
            a.v[4]&=7;
            if (pattern==0) {
                if (i==0) a=zero;
                else if (i<=131) { a=zero;a.v[(i-1)/32]=1u<<((i-1)%32); }
                else if (i==132) a=P131{{~0u,~0u,~0u,~0u,7u}};
            } else if (pattern==1 && (i%3==0 || i%32==0)) a=zero;
            else if (pattern==2) a=zero;
            else if (pattern==3) a=one;
            input[i]=a;
        }
        P131 *di,*dout;checked(cudaMalloc(&di,input.size()*sizeof(P131)));checked(cudaMalloc(&dout,output.size()*sizeof(P131)));
        checked(cudaMemcpy(di,input.data(),input.size()*sizeof(P131),cudaMemcpyHostToDevice));
        checked(cudaMemcpy(dout,output.data(),output.size()*sizeof(P131),cudaMemcpyHostToDevice));
        if (paired) probe<true><<<(n+255)/256? (n+255)/256:1,256>>>(di,dout,n,pitch,guard);
        else probe<false><<<(n+255)/256? (n+255)/256:1,256>>>(di,dout,n,pitch,guard);
        checked(cudaGetLastError());checked(cudaDeviceSynchronize());
        checked(cudaMemcpy(output.data(),dout,output.size()*sizeof(P131),cudaMemcpyDeviceToHost));
        checked(cudaFree(di));checked(cudaFree(dout));
        for (int i=0;i<2*pitch;++i) {
            int local=i%pitch-guard;
            if (local<0 || local>=n) {
                if (!same(output[i],canary)) { fprintf(stderr,"Canary changed: pattern=%d n=%d i=%d\n",pattern,n,i);return 1; }
            } else {
                const P131 a=input[local],got=output[i];
                const int peer = (local & ~255) | ((local & 255) ^ 128);
                const bool iszero=same(a,zero) || (paired && peer<n && same(input[peer],zero));
                const bool valid=(got.v[4]&~7u)==0 && (iszero ? same(got,zero) : i>=pitch ? same(got,a) : same(multiplyReference(a,got),one));
                if (!valid) { fprintf(stderr,"Inverse mismatch: paired=%d pattern=%d n=%d i=%d\n",paired,pattern,n,i);return 1; }
                ++values;if (iszero) ++zeros;
            }
        }
        ++scenarios;
    }
    printf("PASS: %zu block-inverse scenarios, %zu outputs, %zu zero outputs, independent product identities, both zero-propagation modes, repeated shared-tree reuse, ragged CTAs and canaries\n",scenarios,values,zeros);
    return 0;
}
