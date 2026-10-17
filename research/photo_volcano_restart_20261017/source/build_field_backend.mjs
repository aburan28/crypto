// Reproducible source instrumentation only; never edits the frozen backend.
import fs from 'node:fs';
import path from 'node:path';
const root=process.argv[2],out=root+'/field_backend';if(fs.existsSync(out))throw Error('refuse backend overwrite');
fs.cpSync(root+'/source/native_snapshot',out,{recursive:true});
const header=`#ifndef CA_PROBE_H
#define CA_PROBE_H
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct ca_probe_counts {
 uint64_t add[3],sub[3],mul_plain[3],mont_mul[3],mont_square[3],mont_redc[3],inverse[3],euclid_divisions[3],complement[3];
} ca_probe_counts;
#ifdef __cplusplus
extern thread_local ca_probe_counts ca_probe_state;
extern thread_local int ca_probe_enabled;
#else
extern _Thread_local ca_probe_counts ca_probe_state;
extern _Thread_local int ca_probe_enabled;
#endif
static inline unsigned ca_probe_slot(uint64_t m){return m==72057980207847001ULL?0U:m==18014494963342201ULL?1U:2U;}
#define CA_PROBE(field,m) do { if(ca_probe_enabled) ca_probe_state.field[ca_probe_slot(m)]++; } while(0)
static inline uint64_t ca_probe_complement(uint64_t m,uint64_t x){CA_PROBE(complement,m);return m-x;}
CA_API void ca_probe_begin(void);
CA_API ca_probe_counts ca_probe_end(void);
#ifdef __cplusplus
}
#endif
#endif
`;
fs.writeFileSync(out+'/include/cryptanalysis/ca_probe.h',header);
fs.writeFileSync(out+'/src/probe.c','#include "cryptanalysis/ca_types.h"\n#include "cryptanalysis/ca_probe.h"\n#include <string.h>\n_Thread_local ca_probe_counts ca_probe_state;\n_Thread_local int ca_probe_enabled;\nvoid ca_probe_begin(void){memset(&ca_probe_state,0,sizeof ca_probe_state);ca_probe_enabled=1;}\nca_probe_counts ca_probe_end(void){ca_probe_enabled=0;return ca_probe_state;}\n');
const changes=[];
function edit(relative,from,to){const file=path.join(out,relative),text=fs.readFileSync(file,'utf8');if(!text.includes(from))throw Error('instrumentation anchor missing '+relative);const count=text.split(from).length-1;fs.writeFileSync(file,text.replaceAll(from,to));changes.push({file:relative,anchor:from,replacements:count});}
const mod='include/cryptanalysis/ca_modarith.h';edit(mod,'typedef unsigned __int128 ca_u128;','#include "ca_probe.h"\n\ntypedef unsigned __int128 ca_u128;');
for(const [name,args,counter,modulus]of [['ca_addmod','uint64_t a, uint64_t b, uint64_t m','add','m'],['ca_submod','uint64_t a, uint64_t b, uint64_t m','sub','m'],['ca_mulmod','uint64_t a, uint64_t b, uint64_t m','mul_plain','m'],['ca_mont_redc','const ca_mont *m, ca_u128 t','mont_redc','m->p'],['ca_mont_mul','const ca_mont *m, uint64_t a, uint64_t b','mont_mul','m->p'],['ca_mont_sqr','const ca_mont *m, uint64_t a','mont_square','m->p']]){const anchor='static inline uint64_t '+name+'('+args+')\n{';edit(mod,anchor,anchor+'\n    CA_PROBE('+counter+','+modulus+');');}
edit('src/modarith.c','uint64_t ca_invmod(uint64_t a, uint64_t m)\n{','uint64_t ca_invmod(uint64_t a, uint64_t m)\n{\n    CA_PROBE(inverse,m);');
edit('src/modarith.c','    while (newr) {','    while (newr) {\n        CA_PROBE(euclid_divisions,m);');
for(const expression of ['FY(a)','FX(a)','y'])edit('src/group_ec.c','g->p - '+expression,'ca_probe_complement(g->p,'+expression+')');
fs.writeFileSync(out+'/INSTRUMENTATION.json',JSON.stringify({base_source:'source/native_snapshot',arithmetic_expressions_preserved:true,changed_files:changes,threads_supported_for_complete_capture:1,slots:['Fp','scalar_mod_r','other_modulus'],counter_boundary:'Executed modular API calls and explicit EC coordinate complements; not hardware instructions or full state streams'},null,2)+'\n');
console.log(JSON.stringify({out,anchors:changes.length}));
