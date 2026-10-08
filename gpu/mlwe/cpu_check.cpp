#include "mlwe.hpp"
#include <algorithm>
#include <iostream>
#include <random>
#include <vector>

int main() {
    using namespace mlwe;
    for (int k : {2,3,4}) {
        const int batch=3;
        std::vector<int> A(matrix_size(k,batch),0),s(vector_size(k,batch),0),e(s.size(),0);
        // x^255 * x = -1 in this ring, which checks the wrap and its sign.
        A[255]=1; s[1]=1;
        e[0]=2;
        auto t=cpu_reference(A,s,e,k,batch);
        if (t[0]!=1 || t[1]!=0) { std::cerr << "negacyclic wrap failed\n"; return 1; }

        // Large residues exercise sums beyond signed 32-bit range.
        std::fill(A.begin(),A.end(),Q-1);
        std::fill(s.begin(),s.end(),Q-1);
        std::fill(e.begin(),e.end(),Q-1);
        t=cpu_reference(A,s,e,k,batch);
        for (int b=0;b<batch;++b)
            for (int row=0;row<k;++row)
                for (int i=0;i<N;++i)
                    if (t[(b*k+row)*N+i]!=canonical(k*(2*i+2-N)-1)) {
                        std::cerr << "large-residue check failed\n"; return 1;
                    }

        std::mt19937 rng(314159+k); // Reproducible test data; not key-generation randomness.
        std::uniform_int_distribution<int> u(0,Q-1), small(-3,3);
        for (int& a:A) a=u(rng);
        for (int& x:s) x=small(rng);
        for (int& x:e) x=small(rng);
        t=cpu_reference(A,s,e,k,batch);
        // Cross-check selected coefficients through the direct defining equation.
        for (int b=0;b<batch;++b)
            for (int row=0;row<k;++row)
                for (int degree : {0,1,127,255}) {
                    std::int64_t sum=e[(b*k+row)*N+degree];
                    for (int col=0;col<k;++col)
                        for (int i=0;i<N;++i) {
                            const int j=(degree-i+N)%N;
                            const int sign=i<=degree ? 1 : -1;
                            sum+=std::int64_t(sign)*A[(std::size_t(b)*k*k+row*k+col)*N+i]*s[(b*k+col)*N+j];
                        }
                    if (t[(b*k+row)*N+degree]!=canonical(sum)) {
                        std::cerr << "coefficient check failed at k=" << k << "\n";
                        return 1;
                    }
                }
        std::cout << "PASS k=" << k << " batch=" << batch << "\n";
    }
}
