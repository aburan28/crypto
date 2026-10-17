#include "cryptanalysis/ca_modarith.h"
#include <stdio.h>
#define FP 72057980207847001ULL
#define R 18014494963342201ULL
int main(void){ca_mont m;ca_mont_init(&m,FP);uint64_t a=ca_mont_to(&m,2),b=ca_mont_to(&m,3);ca_probe_begin();uint64_t sum=ca_addmod(1,2,FP),diff=ca_submod(8,4,R),prod=ca_mulmod(3,4,97),z=ca_mont_mul(&m,a,b),s=ca_mont_sqr(&m,a),iv=ca_mont_inv(&m,a);ca_probe_counts c=ca_probe_end();if(sum!=3||diff!=4||prod!=12||ca_mont_from(&m,z)!=6||ca_mont_from(&m,s)!=4||ca_mulmod(ca_mont_from(&m,iv),2,FP)!=1)return 1;if(c.add[0]!=1||c.sub[1]!=1||c.mul_plain[2]!=1||c.mont_mul[0]!=2||c.mont_square[0]!=1||c.mont_redc[0]!=3||c.inverse[0]!=1||!c.euclid_divisions[0])return 2;ca_probe_begin();c=ca_probe_end();if(c.mont_mul[0]||c.inverse[0])return 3;puts("probe classification, nesting, arithmetic and reset checks passed");return 0;}
