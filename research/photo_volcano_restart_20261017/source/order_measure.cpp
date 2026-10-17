#include "ec.hpp"
#include <filesystem>
#include <time.h>

struct Counts {
    U field_add_sub=0,field_multiply_including_square=0,field_square=0,field_inverse=0;
    U group_op_calls=0,nontrivial_additions=0,nontrivial_doublings=0,identity_shortcuts=0,inverse_pair_results=0,scalar_calls=0;
    json wire()const{return {{"field_add_sub",field_add_sub},{"field_multiply_including_square",field_multiply_including_square},{"field_square",field_square},{"field_multiply_excluding_square",field_multiply_including_square-field_square},{"field_inverse",field_inverse},{"group_op_calls",group_op_calls},{"nontrivial_additions",nontrivial_additions},{"nontrivial_doublings",nontrivial_doublings},{"identity_shortcuts",identity_shortcuts},{"inverse_pair_results",inverse_pair_results},{"scalar_calls",scalar_calls}};}
};
struct MeasuredEC {
    U a;Counts&c;
    U plus(U x,U y){c.field_add_sub++;return add(x,y);}
    U minus(U x,U y){c.field_add_sub++;return sub(x,y);}
    U times(U x,U y){c.field_multiply_including_square++;return mul(x,y);}
    U square(U x){c.field_square++;return times(x,x);}
    U inverse(U x){need(x!=0,"order inverse zero");c.field_inverse++;U r=1,n=PRIME-2;for(;n;n>>=1,x=square(x))if(n&1)r=times(r,x);return r;}
    Point op(Point p,Point q){
        c.group_op_calls++;
        if(p.infinity||q.infinity){c.identity_shortcuts++;return p.infinity?q:p;}
        if(p.x==q.x&&plus(p.y,q.y)==0){c.inverse_pair_results++;return {};}
        U slope;
        if(p.x==q.x){c.nontrivial_doublings++;slope=times(plus(times(3,square(p.x)),a),inverse(times(2,p.y)));}
        else{c.nontrivial_additions++;slope=times(minus(q.y,p.y),inverse(minus(q.x,p.x)));}
        U x=minus(minus(square(slope),p.x),q.x);
        return {x,minus(times(slope,minus(p.x,x)),p.y),false};
    }
    Point scalar(U k,Point p){c.scalar_calls++;Point result;while(k){if(k&1)result=op(result,p);k>>=1;if(k)p=op(p,p);}return result;}
    U order(Point p){U n=TOTAL;for(U q:{U(2),ORDER})while(n%q==0){if(scalar(n/q,p).infinity)n/=q;else break;}return n;}
};
U process_cpu_ns(){timespec t{};need(clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&t)==0,"CPU clock");return U(t.tv_sec)*1000000000+U(t.tv_nsec);}

int main(int argc,char**argv){try{
    need(argc==3,"ROOT NEW_TAG");std::string root=argv[1],tag=argv[2];
    std::filesystem::path out=root+"/runs/order_"+tag;need(!std::filesystem::exists(out),"refuse overwrite");std::filesystem::create_directories(out);
    auto family=load(root+"/data/family.json"),fixtures=load(root+"/data/fixtures.json");auto binding=load(root+"/ORDER_FROZEN.json");
    need(prime(PRIME)&&prime(ORDER)&&TOTAL==4*ORDER,"known exact order factorization");
    // Independent expected-order controls include the identity, 2-torsion,
    // r-torsion and a point of order 2r. No random tape is used by this phase.
    const auto source=family["nodes"][0];EC reference{val(source["a"]),val(source["b"])};
    Point P=point(fixtures["0"][0]["P"]),T{val(family["edges"][0]["kernel_x"]),0,false};
    std::vector<std::pair<Point,U>>controls={{Point{},1},{T,2},{P,ORDER},{reference.op(P,T),2*ORDER}};
    json control_records=json::array();
    for(const auto&[point,expected]:controls){Counts c;MeasuredEC ec{reference.a,c};need(reference.on(point),"control membership");U result=ec.order(point);need(result==expected,"control order");control_records.push_back({{"expected_order",std::to_string(expected)},{"returned_order",std::to_string(result)},{"counts",c.wire()}});}
    save((out/"controls.json").string(),{{"passed",true},{"controls",control_records}});
    U rows=0;for(const auto&node:family["nodes"]){int cid=node["id"];EC check{val(node["a"]),val(node["b"])};std::ofstream log(out/("curve_"+std::to_string(cid)+".jsonl"));need(bool(log),"order journal");
        for(const auto&fixture:fixtures[std::to_string(cid)])for(const std::string role:{"P","Q"}){
            Point point=::point(fixture[role]);need(!point.infinity&&check.on(point),"fixture membership");Counts counts;MeasuredEC ec{check.a,counts};
            U wall=ns(),cpu=process_cpu_ns();U result=ec.order(point);cpu=process_cpu_ns()-cpu;wall=ns()-wall;
            need(result==ORDER&&check.scalar(result,point).infinity,"independent point order");
            log<<json({{"curve",cid},{"instance",fixture["id"]},{"generator",fixture["generator"]},{"key_index",fixture["key_index"]},{"role",role},{"returned_order",std::to_string(result)},{"wall_ns",wall},{"process_cpu_ns",cpu},{"counters",counts.wire()},{"curve_order_and_factorization_precomputed",true},{"validation_excluded_from_interval",true},{"verified",true}}).dump()<<'\n';rows++;
        }
    }
    save((out/"complete.json").string(),{{"passed",true},{"records",rows},{"expected_records",9200},{"curves",46},{"binding",binding},{"timing_class","CPU-only L0, concurrent BSGS work"},{"counter_boundary","Reference generic point-order algorithm, supplied exact curve order and factorization; input validation and factorization are outside the interval"},{"square_accounting","field_square is a subset of field_multiply_including_square; do not sum both"},{"full_rho_or_bsgs_field_counters",false}});
    std::cout<<rows<<" point-order measurements verified\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
