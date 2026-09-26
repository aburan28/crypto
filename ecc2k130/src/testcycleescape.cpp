// Regression for cycle-entry-dependent exits. No GPU or secret recovery.
#define ECC_NO_CUDA
#define ECC_WALK_TABLE 1
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "../include/curveparams.h"
#include "../include/solver.h"
#include "../include/packedtablewalk.cuh"
using namespace eccPacked131;
static void require(bool ok,const char *why) {
    if(!ok){std::fprintf(stderr,"FAIL: %s\n",why);std::exit(1);}
}
struct GraphOps {
    const unsigned *tags; int n; bool closes; int dp = -1;
    bool distinguished(int p) const {return p==dp;}
    unsigned tag(int p) const {return tags[p%n];}
    bool next(int p,unsigned t,int *q)const {
        require(t==tag(p),"probe uses raw steps");
        *q=closes?(p+1)%n:p+1;return true;
    }
    bool equal(int a,int b)const{return a==b;}
    bool less(int a,int b)const{return a<b;}
};
static void graphCase(const std::vector<unsigned>&tags) {
    GraphOps ops{tags.data(),int(tags.size()),true};
    int anchor=-1;
    for(int entry=0;entry<ops.n;++entry) for(int cold=0;cold<2;++cold) {
        int p=entry; unsigned long long hist=ECC_HIST_EMPTY;
        if(!cold)for(int j=4;j>=1;--j)hist=eccHistPush(hist,tags[(entry+4*ops.n-j)%ops.n]);
        bool escaped=false;
        for(int step=0;step<ops.n+4;++step){
            unsigned raw=ops.tag(p), t=raw;
            if(eccTagFruitless(raw,hist,131))t=eccCycleAnchorTag(p,raw,ops,131,8);
            if(t!=raw){
                if(anchor<0)anchor=p;
                require(anchor==p,"all entry phases choose one exit");escaped=true;break;
            }
            hist=eccHistPush(hist,t);p=(p+1)%ops.n;
        }
        require(escaped,"cold and warm histories escape within one lap plus four steps");
    }
    for (int dp=0;dp<ops.n;++dp) {
        GraphOps reporting{tags.data(),int(tags.size()),true,dp};
        for (int entry=0;entry<ops.n;++entry)
            require(eccCycleAnchorTag(entry,tags[entry],reporting,131,8)==tags[entry],
                    "never skip a distinguished point inside a raw cycle");
    }
    GraphOps open{tags.data(),int(tags.size()),false};
    require(eccCycleAnchorTag(0,tags[0],open,131,8)==tags[0],"false hint never redirects an open path");
}
static P131 pack(const unsigned long long *p){P131 x;for(int i=0;i<5;++i)x.v[i]=uint32_t(p[i/2]>>(32*(i&1)));return x;}
int main(){
    for(int n=2;n<=8;n+=2){
        std::vector<unsigned> tags;
        for(int i=0;i<n/2;++i)tags.push_back(eccTag(i,7*i,0));
        for(int i=0;i<n/2;++i)tags.push_back(eccTag(i,7*i,1));
        graphCase(tags);
    }
    graphCase({eccTag(0,0,0),eccTag(0,0,0),eccTag(0,1,0),eccTag(0,2,0)});
    graphCase({eccTag(0,130,0),eccTag(0,130,0),eccTag(0,0,1),eccTag(0,2,1)});
    using R=Ref<CfgF131>;Solver<CfgF131>s;
    s.setup(eccF131::PX,eccF131::PY,eccF131::QX,eccF131::QY,eccF131::ELL_DEC,eccF131::S_DEC,34,1ull<<32);
    auto p=R::scalarMul(s.basis,u192_from(1184));
    unsigned raw=s.walk.rawTag(p,R::weight(p.x));
    require(raw==0x3e0,"frozen F131 fixture tag");
    auto q=R::addPt(p,s.walk.addend(raw));
    require(s.walk.rawTag(q,R::weight(q.x))==(raw^ECC_TAG_EPS),"frozen raw two-cycle");
    std::vector<uint32_t> table(TW_WORDS);twFillConsts(s.walk,table.data());
    R::Elem exitKey{};bool haveExit=false;int exits=0, comparisons=0;
    for(int phase=0;phase<2;++phase)for(int transform=0;transform<3;++transform){
        auto cur=phase?q:p;
        if(transform==1)cur=R::frob(cur,17);
        if(transform==2)cur=R::neg(cur);
        u64 hist=ECC_HIST_EMPTY;
        bool escaped=false;
        for(int step=0;step<8;++step){
            const int hw=R::weight(cur.x);const unsigned r=s.walk.rawTag(cur,hw);
            const unsigned host=s.walk.resolveTag(cur,r,hist);
            u64 deviceHist=hist;
            const uint32_t *sel=table.data()+TW_SEL0;
            const unsigned device=twSelect(pack(cur.x.v),toPolynomial131(pack(cur.y.v)),hw,&deviceHist,sel,table.data(),34);
            require(host==device,"packed/reference exit tag agreement");
            const auto next=R::addPt(cur,s.walk.addend(host));
            if(host!=r){
                auto key=R::canonical(next.x);
                if(!haveExit){exitKey=key;haveExit=true;}
                require(key==exitKey,"both entry phases and automorphisms use the same exit orbit");
                escaped=true;++exits;
            }
            hist=eccHistPush(hist,host);
            require(hist==deviceHist,"packed/reference history agreement");
            require(R::onCurve(next),"escape stays on curve");cur=next;++comparisons;
            if(escaped)break;
        }
        require(escaped,"actual F131 cycle escaped");
    }
    // An arbitrary inverse hint previously changed every point's branch.
    auto fresh=s.basis;u64 hint=eccHistPush(ECC_HIST_EMPTY,s.walk.rawTag(fresh,R::weight(fresh.x))^ECC_TAG_EPS);
    auto r=s.walk.rawTag(fresh,R::weight(fresh.x));
    require(s.walk.resolveTag(fresh,r,hint)==r,"unrelated history does not force an exit");
    printf("PASS: graph entry phases; F131 exits %d; packed/reference comparisons %d; pivot bytes %d; hybrid %d\n",exits,comparisons,ECC_TABLE_PIVOT_BYTES,ECC_TABLE_ADDEND_GLOBAL);
}
