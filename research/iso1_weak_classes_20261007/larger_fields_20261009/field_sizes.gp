\\ Exact field cardinalities and prime-degree work counts; logarithms are display data.
C(p,n)={my(N=(p^(2*n)-1)/(p^2-1),A=(p^n-1)/(p-1),B=(p^n+1)/(p+1),e=(p^2%n==1));(N+A+B-3+2*e*(n-1)^2)/(4*n)};
cases=[[13,3],[257,3],[509,3],[1009,3],[2003,3],[8191,3],[65537,3],[1048583,3],[4294967311,3],[4398046511119,3],[13,5],[257,5],[1009,5],[65537,5],[13,7],[257,7],[1009,7],[65537,7]];
for(i=1,#cases,p=cases[i][1];n=cases[i][2];Q=p^(2*n);print(p,"|",n,"|",Q,"|",strprintf("%.6f",log(Q)/log(2)),"|",C(p,n),"|",strprintf("%.6f",log(C(p,n))/log(10))));
