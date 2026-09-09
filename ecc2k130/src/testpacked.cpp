#include "../include/curveparams.h"
#include "../include/packed131.h"
#include <cstdio>
using R = Ref<CfgF131>;
using P = eccPacked131::P131;
static unsigned long long state = 0x131ab123456789ULL;
static unsigned long long randomWord() { state ^= state << 13; state ^= state >> 7; state ^= state << 17; return state; }
static P pack(R::Elem a) { P p; for (int i=0;i<5;i++) p.v[i]=unsigned(a.v[i/2]>>(32*(i&1))); return p; }
static R::Elem unpack(P p) { unsigned long long a[3]={p.v[0]|(static_cast<unsigned long long>(p.v[1])<<32),p.v[2]|(static_cast<unsigned long long>(p.v[3])<<32),p.v[4]};return R::fromLimbs(a); }
int main() {
    // Check every pair of unit coefficients, including the top-three-bit
    // product and the output fold. Dense cases below also exercise carries in
    // the integer-mask carryless primitive.
    for (int i=0;i<131;i++) for (int j=0;j<131;j++) {
        unsigned long long av[3]={},bv[3]={};
        av[i/64]=1ull<<(i%64);bv[j/64]=1ull<<(j%64);
        auto a=R::fromLimbs(av),b=R::fromLimbs(bv);
        auto pa=eccPacked131::toPolynomial131(pack(a)),pb=eccPacked131::toPolynomial131(pack(b));
        auto product=eccPacked131::fromPolynomial131(eccPacked131::mulPolynomial131(pa,pb));
        if (unpack(eccPacked131::mul131(pack(a),pack(b)))!=R::mul(a,b) || unpack(product)!=R::mul(a,b)) {
            printf("packed basis-pair mismatch at %d,%d\n",i,j);return 1;
        }
    }
    for (int test=0;test<160;test++) {
        unsigned long long av[3]={randomWord(),randomWord(),randomWord()&7},bv[3]={randomWord(),randomWord(),randomWord()&7};
        if (test<131) {av[0]=av[1]=av[2]=0;av[test/64]=1ull<<(test%64);}
        auto a=R::fromLimbs(av),b=R::fromLimbs(bv);auto pa=pack(a),pb=pack(b);
        auto ap=eccPacked131::toPolynomial131(pa),bp=eccPacked131::toPolynomial131(pb);
        if (unpack(eccPacked131::fromPolynomial131(ap))!=a ||
            unpack(eccPacked131::fromPolynomial131(eccPacked131::mulPolynomial131(ap,bp)))!=R::mul(a,b)) {
            printf("polynomial field mismatch at case %d\n",test);return 1;
        }
        if (unpack(eccPacked131::mul131(pa,pb))!=R::mul(a,b) ||
            unpack(eccPacked131::sqr131(pa))!=R::sqr(a) ||
            unpack(eccPacked131::inv131(pa))!=R::inv(a)) {
            printf("packed field mismatch at case %d\n",test);return 1;
        }
        const int powers[]={0,1,2,3,4,5,6,7,8,9,10,16,32,65,130,131};
        for (int j:powers) if (unpack(eccPacked131::sigma131(pa,j))!=R::sigma(a,j)) {
            printf("packed Frobenius mismatch at case %d, power %d\n",test,j);return 1;
        }
    }
    P zero{};auto one=pack(R::one());
    if (unpack(eccPacked131::mul131(zero,one))!=R::zero() ||
        unpack(eccPacked131::mul131(one,one))!=R::one()) return 1;
    puts("PASS: packed multiplication, squaring, inversion, walk and inversion Frobenius powers against independent reference");
}
