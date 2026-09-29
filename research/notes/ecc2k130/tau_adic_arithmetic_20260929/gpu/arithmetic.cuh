// Portable known-scalar reference arithmetic. No walk or discrete-log solver.
// FIELD_M is 83 or 131; FAST_SQUARE selects the independently checked variant.
typedef unsigned long long U;
#ifndef FIELD_M
#error FIELD_M must be defined
#endif
#ifndef FAST_SQUARE
#define FAST_SQUARE 0
#endif
#if FIELD_M != 83 && FIELD_M != 131
#error Unsupported fixture field
#endif
#ifdef __CUDACC__
#define HD __host__ __device__ __forceinline__
#else
#define HD inline
#endif
static const int W = (FIELD_M+63)/64;
static const U TOP = (1ULL << (FIELD_M%64))-1;
static const U LOW = FIELD_M == 83 ? ((1ULL<<45)|7) : ((1ULL<<13)|7);
struct F { U v[W]; };
struct P { F x, y; bool inf; };
HD F zero() { F a; for(int j=0;j<W;++j) a.v[j]=0; return a; }
HD bool iszero(F a) { U x=0; for(int j=0;j<W;++j)x|=a.v[j]; return x==0; }
HD bool equal(F a,F b) { U x=0; for(int j=0;j<W;++j)x|=a.v[j]^b.v[j];return x==0; }
HD F plus(F a,F b) { for(int j=0;j<W;++j)a.v[j]^=b.v[j];return a; }
HD F oneplus(F a) { a.v[0]^=1;return a; }
HD F times(F a,F b) {
    F out=zero();
    for(int i=0;i<FIELD_M;++i) {
        const U mask=0ULL-((b.v[i/64]>>(i%64))&1ULL);
        for(int j=0;j<W;++j)out.v[j]^=a.v[j]&mask;
        const U carry=(a.v[W-1]>>((FIELD_M-1)%64))&1;
        for(int j=W-1;j>0;--j)a.v[j]=(a.v[j]<<1)|(a.v[j-1]>>63);
        a.v[0]<<=1; a.v[W-1]&=TOP; a.v[0]^=LOW*(U)carry;
    }
    return out;
}
HD F square(F a,const U* columns) {
#if FAST_SQUARE
    F out=zero();
    for(int i=0;i<FIELD_M;++i) {
        const U mask=0ULL-((a.v[i/64]>>(i%64))&1ULL);
        for(int j=0;j<W;++j)out.v[j]^=columns[i*W+j]&mask;
    }
    return out;
#else
    (void)columns; return times(a,a);
#endif
}
HD F inverse(F a,const U* columns) {
    // Itoh-Tsujii: construct a^(2^(m-1)-1), then square.
    F t=a; int k=1; const int e=FIELD_M-1;
    int top=0; while((1<<(top+1))<=e)++top;
    for(int bit=top-1;bit>=0;--bit) {
        F u=t; for(int j=0;j<k;++j)u=square(u,columns);
        t=times(u,t); k*=2;
        if((e>>bit)&1) { t=times(square(t,columns),a); ++k; }
    }
    return square(t,columns);
}
HD P infinity() { P p; p.x=zero();p.y=zero();p.inf=true;return p; }
HD P negate(P p) { if(!p.inf)p.y=plus(p.x,p.y);return p; }
HD P twice(P p,const U* columns) {
    if(p.inf||iszero(p.x))return infinity();
    F lambda=plus(p.x,times(p.y,inverse(p.x,columns)));
    P r; r.inf=false;r.x=plus(square(lambda,columns),lambda);
    r.y=plus(square(p.x,columns),times(oneplus(lambda),r.x));return r;
}
HD P add(P p,P q,const U* columns) {
    if(p.inf)return q;if(q.inf)return p;
    if(equal(p.x,q.x))return equal(p.y,q.y)?twice(p,columns):infinity();
    F d=plus(p.x,q.x),lambda=times(plus(p.y,q.y),inverse(d,columns));
    P r; r.inf=false;r.x=plus(plus(square(lambda,columns),lambda),d);
    r.y=plus(plus(times(lambda,plus(p.x,r.x)),r.x),p.y);return r;
}
HD P frob(P p,const U* columns) {
    if(!p.inf){p.x=square(p.x,columns);p.y=square(p.y,columns);}return p;
}
HD void evaluate_one(const U* points,const signed char* digits,const int* lengths,
                     int n,const U* columns,U* output,int use_tau,int i) {
    P p;for(int j=0;j<W;++j){p.x.v[j]=points[(2*j)*n+i];p.y.v[j]=points[(2*j+1)*n+i];}
    p.inf=points[2*W*n+i]!=0;
    P r=infinity(),neg=negate(p);
    for(int bit=lengths[i]-1;bit>=0;--bit) {
        r=use_tau?frob(r,columns):twice(r,columns);
        const int d=digits[bit*n+i];if(d)r=add(r,d==1?p:neg,columns);
    }
    for(int j=0;j<W;++j){output[(2*j)*n+i]=r.inf?0:r.x.v[j];output[(2*j+1)*n+i]=r.inf?0:r.y.v[j];}
    output[2*W*n+i]=r.inf?1:0;
}
#ifdef __CUDACC__
extern "C" __global__ void evaluate(const U* points,const signed char* digits,
    const int* lengths,int n,const U* columns,U* output,int use_tau) {
    const int i=blockIdx.x*blockDim.x+threadIdx.x;
    if(i<n)evaluate_one(points,digits,lengths,n,columns,output,use_tau,i);
}
#else
extern "C" void host_evaluate(const U* points,const signed char* digits,
    const int* lengths,int n,const U* columns,U* output,int use_tau) {
    for(int i=0;i<n;++i)evaluate_one(points,digits,lengths,n,columns,output,use_tau,i);
}
extern "C" void host_field(const U* input,const U* columns,U* sq,U* inv) {
    F p;for(int j=0;j<W;++j)p.v[j]=input[j];
    F a=square(p,columns),b=inverse(p,columns);
    for(int j=0;j<W;++j){sq[j]=a.v[j];inv[j]=b.v[j];}
}
#endif
