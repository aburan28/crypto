\\ Attempt the existing exact CM support path on one admitted independent source.
\\ Parent driver enforces a 60-second total cap and retains the last printed stage.
{
  my(p=4294967291,t=26895410567024103251795854646,Q=p^6,delta=t^2-4*Q);
  print("CM_PREFLIGHT|p=",p,"|trace=",t,"|delta=",delta,"|stage=FUNDAMENTAL_DISCRIMINANT");
  my(dk=coredisc(delta),f=sqrtint(delta/dk));
  if(dk*f^2!=delta||f%4,error("CM support preflight conductor check failed"));
  print("CM_PREFLIGHT|D_K=",dk,"|f_pi=",f,"|stage=MAXIMAL_ORDER_CLASS_POLYNOMIAL");
  my(H=polclass(dk));
  print("CM_PREFLIGHT|degree=",poldegree(H),"|stage=FIRST_ORDER_COMPLETE");
}
quit();
