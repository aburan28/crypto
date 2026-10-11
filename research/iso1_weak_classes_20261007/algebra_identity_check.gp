\\ Exact expansion over Z, independent of native modular evaluations.
quotient = x^2 + 2*(1+s^2)*x + (1-s^2)^2 - (x+(1+s)^2)*(x+(1-s)^2);
if(quotient != 0, error("quotient factorization failed"));
U = x^2 + a*x + b;
isogeny = U^2 - 2*a*x*U + (a^2-4*b)*x^2 - (x^2-b)^2;
if(isogeny != 0, error("cleared isogeny equation failed"));
print("integer-polynomial quotient factorization: exact zero");
print("integer-polynomial cleared 2-isogeny equation: exact zero");
