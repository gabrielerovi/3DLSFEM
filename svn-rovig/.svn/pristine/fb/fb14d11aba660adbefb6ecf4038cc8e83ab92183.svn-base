function x=solve_2x2system(A,b)

detA=A(1,1)*A(2,2)-A(1,2)*A(2,1);
A_minus1=1.0/detA *[A(2,2) -A(1,2);-A(2,1) A(1,1)];

x=A_minus1*b;

end