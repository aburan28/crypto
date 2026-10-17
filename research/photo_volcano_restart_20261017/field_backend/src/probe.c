#include "cryptanalysis/ca_types.h"
#include "cryptanalysis/ca_probe.h"
#include <string.h>
_Thread_local ca_probe_counts ca_probe_state;
_Thread_local int ca_probe_enabled;
void ca_probe_begin(void){memset(&ca_probe_state,0,sizeof ca_probe_state);ca_probe_enabled=1;}
ca_probe_counts ca_probe_end(void){ca_probe_enabled=0;return ca_probe_state;}
