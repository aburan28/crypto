// EXPERIMENTAL MODIFIED CMS FUNCTIONAL PRIMARY-ONLY CDCL FORK.
// Linked standalone to a modified public source snapshot, not the pinned wheel.
#include "cryptominisat.h"
#include <atomic>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
thread_local std::string error;
struct Context {
    std::atomic<bool> interrupted{false};
    std::unique_ptr<CMSat::SATSolver> solver;
    uint32_t variables;
    explicit Context(uint32_t n, uint32_t primary) : variables(n) {
        solver.reset(new CMSat::SATSolver(nullptr, &interrupted));
        solver->set_num_threads(1);
        solver->set_seed(0);
        solver->set_no_simplify();
        solver->set_no_simplify_at_startup();
        solver->set_no_bve();
        solver->set_no_bva();
        solver->set_no_equivalent_lit_replacement();
        solver->set_renumber(false);
        solver->set_sls(0);
        solver->set_find_xors(false);
        solver->set_xor_detach(false);
        solver->new_vars(n);
        // new_vars is queued; implied_by does not flush it. This public query
        // explicitly flushes pending variable/clause additions in one-thread mode.
        solver->get_lit_incidence();
        solver->set_functional_primary_mask(primary);
    }
};
CMSat::Lit literal(int32_t value, uint32_t n) {
    const int64_t absolute = value < 0 ? -int64_t(value) : value;
    if (absolute < 1 || absolute > n) throw std::invalid_argument("literal out of declared range");
    return CMSat::Lit(uint32_t(absolute-1), value < 0);
}
std::vector<CMSat::Lit> literals(const int32_t* values, size_t count, uint32_t n) {
    std::vector<CMSat::Lit> result;
    result.reserve(count);
    for (size_t i=0; i<count; ++i) result.push_back(literal(values[i],n));
    return result;
}
int32_t encoded(CMSat::Lit value) {
    return value.sign() ? -int32_t(value.var()+1) : int32_t(value.var()+1);
}
}
extern "C" {
const char* pd_error() { return error.c_str(); }
const char* pd_version() { return CMSat::SATSolver::get_version(); }
const char* pd_git_sha() { return CMSat::SATSolver::get_version_sha1(); }
void* pd_create(uint32_t variables,uint32_t primary) {
    try { error.clear(); return new Context(variables,primary); }
    catch (const std::exception& e) { error=e.what(); return nullptr; }
}
void pd_destroy(void* pointer) { delete static_cast<Context*>(pointer); }
void pd_cancel(void* pointer) { static_cast<Context*>(pointer)->interrupted.store(true); }
int pd_add_clause(void* pointer,const int32_t* values,size_t count) {
    try { auto& c=*static_cast<Context*>(pointer); return c.solver->add_clause(literals(values,count,c.variables)) ? 1:0; }
    catch (const std::exception& e) { error=e.what(); return -2; }
}
int pd_add_xor(void* pointer,const uint32_t* values,size_t count,int parity) {
    try {
        auto& c=*static_cast<Context*>(pointer);
        if (count>16384) throw std::invalid_argument("declared long XOR row cap exceeded");
        std::vector<unsigned> row;
        for(size_t i=0;i<count;++i) {
            if(values[i]<1 || values[i]>c.variables) throw std::invalid_argument("XOR variable out of range");
            row.push_back(values[i]-1);
        }
        return c.solver->add_xor_clause(row,bool(parity)) ? 1:0;
    } catch (const std::exception& e) { error=e.what(); return -2; }
}
int pd_implied(void*,const int32_t*,size_t,int32_t*,size_t,size_t*) {
 error="implied_by is disabled in the long-XOR interface because native cut variables are not declared coordinates";
 return -2;
}
int pd_solve(void* pointer,const int32_t* prefix,size_t count,double seconds,uint64_t conflicts,int8_t* model,size_t capacity) {
    try {
        auto& c=*static_cast<Context*>(pointer);
        auto assumptions=literals(prefix,count,c.variables);
        c.solver->set_max_time(seconds);
        c.solver->set_max_confl(conflicts);
        auto verdict=c.solver->solve(&assumptions);
        if(verdict==CMSat::l_Undef) return -1;
        if(verdict==CMSat::l_False) return 0;
        const auto& state=c.solver->get_model();
        if(state.size()<c.variables || capacity<c.variables+1) throw std::length_error("full native model is missing declared coordinates");
        model[0]=-1;
        for(uint32_t i=0;i<c.variables;++i) {
            if(state[i]==CMSat::l_Undef) throw std::runtime_error("native leaf model has an undefined functional coordinate");
            model[i+1]=state[i]==CMSat::l_True ? 1:0;
        }
        return 1;
    } catch (const std::exception& e) { error=e.what(); return -2; }
}
int pd_flush(void* pointer) {
 try { static_cast<Context*>(pointer)->solver->get_lit_incidence(); return 1; }
 catch (const std::exception& e) { error=e.what(); return -2; }
}
void pd_stats(void* pointer,uint64_t* output) {
    auto& s=*static_cast<Context*>(pointer)->solver;
    output[0]=s.get_sum_conflicts();output[1]=s.get_sum_decisions();output[2]=s.get_sum_propagations();
    s.get_functional_primary_stats(output+3);
}
}
