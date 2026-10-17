#define FINAL_AGGREGATE_LIBRARY
#include "final_aggregate.cpp"
int main(){
 completion_gate(55,50,92,30,5);
 size_t rejected=0;
 for(auto c:std::vector<std::array<size_t,5>>{{55,53,116,30,5},{55,50,92,29,5},{54,50,92,30,5},{55,50,92,30,4}}){try{completion_gate(c[0],c[1],c[2],c[3],c[4]);}catch(...){rejected++;}}
 require(rejected==4,"Coverage guard");
 J a=Parser{"{\"nested\":{\"work\":1}}"}.value(),b=Parser{"{\"nested\":{\"work\":2}}"}.value();
 require(!same(a,b)&&same(a,a),"Recursive A/A");
 auto p=files::temp_directory_path()/("n83-binding-test-"+std::to_string(getpid()));
 create_new(p,"original\n");J hash;hash.type="string";hash.scalar=digest(p);sealed(p,hash);
 {std::ofstream f(p);f<<"changed\n";}
 bool tampered=false;try{sealed(p,hash);}catch(...){tampered=true;}
 require(tampered,"Modified raw evidence accepted");files::remove(p);
 std::cout<<"PASS actual/planned coverage, recursive A/A tamper and raw SHA binding guards\n";
}
