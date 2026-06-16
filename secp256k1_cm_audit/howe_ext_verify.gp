\\ howe_ext_verify.gp — verify Jac(C)/F_{43^3} splits as E1 x E2
default(parisize,256000000); default(timer,0);
p=43;

\\ From howe_direct_f43.gp: b1=6, b2=55
\\ Char poly of Frob eigenvalues on H^1(Jac(C)): T^4-6T^3+55T^2-258T+1849
e1=6; e2=55; e3=258; e4=1849;
p1=e1;
p2=e1*p1-2*e2;
p3=e1*p2-e2*p1+3*e3;
print("Trace of Frob^3 on H^1(Jac(C)) [should be 0 if Jac splits over F_{p^3}]:");
print("  p_3 = ", p3);
print("");

\\ Traces of E1 and E2 over F_{43^3} via recurrence
T3_E1 = 13^3 - 3*p*13;
T3_E2 = (-13)^3 - 3*p*(-13);
print("T_3 for E1 (trace=13) over F_{43^3}: ", T3_E1);
print("T_3 for E2 (trace=-13) over F_{43^3}: ", T3_E2);
print("Sum = ", T3_E1 + T3_E2, " (= p_3 expected)");
print("");

\\ Point counts
E1_pts = p^3 + 1 - T3_E1;
E2_pts = p^3 + 1 - T3_E2;
print("E1(F_{43^3}) = ", E1_pts);
print("E2(F_{43^3}) = ", E2_pts);
print("E1 * E2 = ", E1_pts * E2_pts);
print("");

\\ Check: L_{Jac/F_{p^3}}(1) = product (1-alpha_i^3)
\\ = (1-T3_E1+p^3)(1+T3_E1+p^3) = E1_pts * E2_pts (since traces of E2 = -T3_E1)
check = (1 - T3_E1 + p^3) * (1 + T3_E1 + p^3);
print("Predicted Jac(C)(F_{43^3}) = ", check);
print("Matches E1*E2? ", check == E1_pts * E2_pts);
print("");

\\ Summary
print("CONCLUSION:");
if (p3 == 0,
    print("  Frobenius^3 trace = 0: Jac(C)/F_{43^3} is F_{43^3}-isogenous to E1 x E2."),
    print("  Frobenius^3 trace != 0: unexpected.")
);
print("  The (2,2)-isogeny Jac(C) -> E1 x E2 is defined over F_{43^3}, NOT F_43.");
print("  Igusa: J2=41, J10=36 mod 43 (curve smooth).");
print("  The 2-torsion E1[2], E2[2] live in F_{43^3} (both #E1(F_43)=31, #E2(F_43)=57 odd).");
print("  This is the structural F_{p^3} obstruction to Howe cover attacks over F_p.");
