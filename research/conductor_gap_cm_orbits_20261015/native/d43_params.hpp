#include <map>
thread_local U CM_P,CM_R,CM_PAIR,CM_LAMBDA,CM_INV_C,CM_INV_CU;
thread_local U CM_X[5],CM_T[5],CM_Y2[5];
struct CMParameters{U p,r,pair,lambda,inv_c,inv_cu,x[5],t[5],y2[5];};
std::map<U,CMParameters> all_cm_parameters;
