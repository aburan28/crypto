\\ Exact integer polynomial expansion for the new theorem records.
halving=(a*b*(a+b))^2-a*b*(a*b+a^2)*(a*b+b^2);
if(halving!=0,error("halving point identity failed"));
orbit=(p^4-p^2)+4*(p^2-1)+12-(p^4+3*p^2+8);
if(orbit!=0,error("orbit numerator identity failed"));
print("rational halving point: exact integer-polynomial zero");
print("absolute-orbit numerator: exact integer-polynomial zero");
