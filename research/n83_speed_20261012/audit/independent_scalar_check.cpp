#include <gmpxx.h>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

using U = unsigned __int128;
using cpp_int = mpz_class;
static constexpr uint64_t R = 8569786107849059ULL;
static uint64_t state = 202610126281ULL;
uint64_t next() {
    state += 0x9e3779b97f4a7c15ULL;
    auto z = state;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
cpp_int big(U x) { return (cpp_int(uint64_t(x >> 64)) << 64) + uint64_t(x); }
U high(U a, U b) {
    U a0=uint64_t(a), a1=a>>64, b0=uint64_t(b), b1=b>>64;
    U p00=a0*b0, p01=a0*b1, p10=a1*b0, p11=a1*b1;
    U carry=(p00>>64)+uint64_t(p01)+uint64_t(p10);
    return p11+(p01>>64)+(p10>>64)+(carry>>64);
}
void require(bool ok, const char *message) { if (!ok) throw std::runtime_error(message); }
int main() {
    try {
        const U max=~U(0);
        const cpp_int B=cpp_int(1)<<128;
        uint64_t products=0, highs=0, adds=0, corrections=0;
        const std::vector<U> words={0,1,U(1)<<64,(U(1)<<64)-1,U(1)<<127,max-1,max};
        auto check_high=[&](U a,U b) {
            require(big(high(a,b))==((big(a)*big(b))>>128),"high128 mismatch"); ++highs;
        };
        for (auto a:words) for(auto b:words) check_high(a,b);
        for(unsigned i=0;i<100000;++i) {
            auto ah=next(), al=next(), bh=next(), bl=next();
            check_high((U(ah)<<64)|al,(U(bh)<<64)|bl);
        }
        std::vector<uint64_t> moduli={2,3,4,5,7,8,15,16,17,65535,65536,65537,R-1,R,R+1,
            (uint64_t(1)<<50)-1,uint64_t(1)<<50,(uint64_t(1)<<50)+1,
            (uint64_t(1)<<63)-1,uint64_t(1)<<63,(uint64_t(1)<<63)+1,UINT64_MAX-1,UINT64_MAX};
        for(unsigned i=0;i<1000;++i) { auto m=next(); if(m<2)m=2; moduli.push_back(m); }
        for(auto m:moduli) {
            const U mu=max/m+U(max%m==m-1);
            require(big(mu)==B/m,"reciprocal mismatch");
            std::vector<uint64_t> residues={0,1,m/2,m-2,m-1};
            auto check_pair=[&](uint64_t a,uint64_t b) {
                const U x=U(a)*b,q=high(x,mu),exact=x/m;
                require(q<=exact && exact-q<=1,"quotient bound mismatch");
                const U residual=x-q*m;
                require(residual<U(m)*2,"two-modulus residual bound mismatch");
                const U reduced=residual>=m?residual-m:residual;
                require(reduced==x%m,"Barrett remainder mismatch"); ++products;
                corrections+=residual>=m;
                auto sum=uint64_t(a+b);
                bool carry=sum<a;
                auto add=(carry||sum>=m)?uint64_t(sum-m):sum;
                require(U(add)==(U(a)+b)%m,"overflow-aware add mismatch"); ++adds;
            };
            for(auto a:residues) for(auto b:residues) check_pair(a,b);
            for(unsigned i=0;i<1000;++i) {
                auto a=next()%m, b=next()%m;
                check_pair(a,b);
            }
        }
        std::cout<<"{\"schema\":\"n83-independent-scalar-audit/v1\",\"status\":\"PASS\","
            "\"seed\":202610126281,\"moduli\":"<<moduli.size()<<",\"high128_checks\":"<<highs
            <<",\"multiply_checks\":"<<products<<",\"add_checks\":"<<adds
            <<",\"one_corrections_observed\":"<<corrections
            <<",\"oracle\":\"GMP arbitrary-precision high product and reciprocal; exact unsigned128 remainder\","
            "\"scope\":\"finite boundary and deterministic random tests, not universal proof or solver trace\"}\n";
        return 0;
    } catch(const std::exception &e) { std::cerr<<e.what()<<"\n"; return 1; }
}
