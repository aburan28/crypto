/* Explicit product-row carry chains for sm_120 and other NVIDIA GPUs.
 * Include this header (rather than nist_fields.cuh directly) in CUDA code.
 */
#ifndef GPU_NIST_PTX_CUH
#define GPU_NIST_PTX_CUH
#include "nist_fields.cuh"
#if defined(__CUDACC__) && NIST_PTX

template<> NF_HD void n_mac_row<8>(uint32_t*t,const uint32_t*a,uint32_t b){
#if defined(__CUDA_ARCH__)
 asm("mad.lo.cc.u32 %0,%10,%18,%0;\n\t"
     "madc.lo.cc.u32 %1,%11,%18,%1;\n\t"
     "madc.lo.cc.u32 %2,%12,%18,%2;\n\t"
     "madc.lo.cc.u32 %3,%13,%18,%3;\n\t"
     "madc.lo.cc.u32 %4,%14,%18,%4;\n\t"
     "madc.lo.cc.u32 %5,%15,%18,%5;\n\t"
     "madc.lo.cc.u32 %6,%16,%18,%6;\n\t"
     "madc.lo.cc.u32 %7,%17,%18,%7;\n\t"
     "addc.cc.u32 %8,%8,0;\n\taddc.u32 %9,%9,0;"
     :"+r"(t[0]),"+r"(t[1]),"+r"(t[2]),"+r"(t[3]),"+r"(t[4]),"+r"(t[5]),"+r"(t[6]),"+r"(t[7]),"+r"(t[8]),"+r"(t[9])
     :"r"(a[0]),"r"(a[1]),"r"(a[2]),"r"(a[3]),"r"(a[4]),"r"(a[5]),"r"(a[6]),"r"(a[7]),"r"(b));
 asm("mad.hi.cc.u32 %1,%10,%18,%1;\n\t"
     "madc.hi.cc.u32 %2,%11,%18,%2;\n\t"
     "madc.hi.cc.u32 %3,%12,%18,%3;\n\t"
     "madc.hi.cc.u32 %4,%13,%18,%4;\n\t"
     "madc.hi.cc.u32 %5,%14,%18,%5;\n\t"
     "madc.hi.cc.u32 %6,%15,%18,%6;\n\t"
     "madc.hi.cc.u32 %7,%16,%18,%7;\n\t"
     "madc.hi.cc.u32 %8,%17,%18,%8;\n\taddc.u32 %9,%9,0;"
     :"+r"(t[0]),"+r"(t[1]),"+r"(t[2]),"+r"(t[3]),"+r"(t[4]),"+r"(t[5]),"+r"(t[6]),"+r"(t[7]),"+r"(t[8]),"+r"(t[9])
     :"r"(a[0]),"r"(a[1]),"r"(a[2]),"r"(a[3]),"r"(a[4]),"r"(a[5]),"r"(a[6]),"r"(a[7]),"r"(b));
#else
 uint64_t c=0;for(int j=0;j<8;j++){c+=(uint64_t)t[j]+(uint64_t)a[j]*b;t[j]=(uint32_t)c;c>>=32;}
 c+=t[8];t[8]=(uint32_t)c;t[9]+=(uint32_t)(c>>32);
#endif
}

template<> NF_HD void n_mac_row<12>(uint32_t*t,const uint32_t*a,uint32_t b){
#if defined(__CUDA_ARCH__)
 asm("mad.lo.cc.u32 %0,%14,%26,%0;\n\t"
     "madc.lo.cc.u32 %1,%15,%26,%1;\n\t"
     "madc.lo.cc.u32 %2,%16,%26,%2;\n\t"
     "madc.lo.cc.u32 %3,%17,%26,%3;\n\t"
     "madc.lo.cc.u32 %4,%18,%26,%4;\n\t"
     "madc.lo.cc.u32 %5,%19,%26,%5;\n\t"
     "madc.lo.cc.u32 %6,%20,%26,%6;\n\t"
     "madc.lo.cc.u32 %7,%21,%26,%7;\n\t"
     "madc.lo.cc.u32 %8,%22,%26,%8;\n\t"
     "madc.lo.cc.u32 %9,%23,%26,%9;\n\t"
     "madc.lo.cc.u32 %10,%24,%26,%10;\n\t"
     "madc.lo.cc.u32 %11,%25,%26,%11;\n\t"
     "addc.cc.u32 %12,%12,0;\n\taddc.u32 %13,%13,0;"
     :"+r"(t[0]),"+r"(t[1]),"+r"(t[2]),"+r"(t[3]),"+r"(t[4]),"+r"(t[5]),"+r"(t[6]),"+r"(t[7]),"+r"(t[8]),"+r"(t[9]),"+r"(t[10]),"+r"(t[11]),"+r"(t[12]),"+r"(t[13])
     :"r"(a[0]),"r"(a[1]),"r"(a[2]),"r"(a[3]),"r"(a[4]),"r"(a[5]),"r"(a[6]),"r"(a[7]),"r"(a[8]),"r"(a[9]),"r"(a[10]),"r"(a[11]),"r"(b));
 asm("mad.hi.cc.u32 %1,%14,%26,%1;\n\t"
     "madc.hi.cc.u32 %2,%15,%26,%2;\n\t"
     "madc.hi.cc.u32 %3,%16,%26,%3;\n\t"
     "madc.hi.cc.u32 %4,%17,%26,%4;\n\t"
     "madc.hi.cc.u32 %5,%18,%26,%5;\n\t"
     "madc.hi.cc.u32 %6,%19,%26,%6;\n\t"
     "madc.hi.cc.u32 %7,%20,%26,%7;\n\t"
     "madc.hi.cc.u32 %8,%21,%26,%8;\n\t"
     "madc.hi.cc.u32 %9,%22,%26,%9;\n\t"
     "madc.hi.cc.u32 %10,%23,%26,%10;\n\t"
     "madc.hi.cc.u32 %11,%24,%26,%11;\n\t"
     "madc.hi.cc.u32 %12,%25,%26,%12;\n\taddc.u32 %13,%13,0;"
     :"+r"(t[0]),"+r"(t[1]),"+r"(t[2]),"+r"(t[3]),"+r"(t[4]),"+r"(t[5]),"+r"(t[6]),"+r"(t[7]),"+r"(t[8]),"+r"(t[9]),"+r"(t[10]),"+r"(t[11]),"+r"(t[12]),"+r"(t[13])
     :"r"(a[0]),"r"(a[1]),"r"(a[2]),"r"(a[3]),"r"(a[4]),"r"(a[5]),"r"(a[6]),"r"(a[7]),"r"(a[8]),"r"(a[9]),"r"(a[10]),"r"(a[11]),"r"(b));
#else
 uint64_t c=0;for(int j=0;j<12;j++){c+=(uint64_t)t[j]+(uint64_t)a[j]*b;t[j]=(uint32_t)c;c>>=32;}
 c+=t[12];t[12]=(uint32_t)c;t[13]+=(uint32_t)(c>>32);
#endif
}
#endif
#endif
