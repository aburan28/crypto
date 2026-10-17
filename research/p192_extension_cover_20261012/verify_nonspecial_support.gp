\\ Independently check the degree-two y-interpolation coefficients in F_p.
read("research/p192_extension_cover_20261012/type_i_nonspecial_support_input.gp");
verify()={
  my(nonzero=0);
  if(P != 6277101735386680763835789423207666416083908700390324961279,
    error("wrong field prime"));
  if(#Y != 160, error("wrong class count"));
  for(i=1,#Y,
    if(poldegree(Y[i]) != 2, error(Str("collinear support at index ",i)));
    if(Mod(polcoeff(Y[i],2),P) == 0, error("zero quadratic coefficient"));
    nonzero++;
  );
  print("NONSPECIAL_SUPPORT_GP_PASS classes=",#Y," nonzero_quadratic=",nonzero);
};
verify();
quit();
