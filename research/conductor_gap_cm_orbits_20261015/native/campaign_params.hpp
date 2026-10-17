#include <map>
thread_local U CM_P,CM_R,CM_PAIR,CM_LAMBDA,CM_INV_C,CM_INV_CU;
thread_local U CM_X[2],CM_T[2],CM_Y2[2];
struct CMParameters{U p,r,pair,lambda,inv_c,inv_cu,x[2],t[2],y2[2];};
std::map<U,CMParameters> all_cm_parameters;
