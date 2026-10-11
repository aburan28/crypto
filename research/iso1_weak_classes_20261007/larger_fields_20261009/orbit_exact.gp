\\ Cleared prime-degree orbit sum, independently expanded over the integers.
check_orbit_identity()={
  my(lhs=2*n*e*(n-1)+2*(A+B-2-e*(n-1))+(N+1-A-B));
  my(rhs=N+A+B-3+2*e*(n-1)^2);
  if(lhs-rhs,error("cleared orbit identity failed"));
  print("prime-degree cleared orbit sum: exact integer-polynomial zero");
};
check_orbit_identity();
