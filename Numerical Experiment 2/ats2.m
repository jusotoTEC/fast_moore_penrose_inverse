function X=ats2(A)
% Method to compute Pseudoinverse of A
% Sintaxis:  X = ats2(A)
%                Input:  Matrix X of size m x n (m>=n)
%                Output: Pseudoinverse X of size n x m
% Reference: Stanimirovic, P. S., Pappas, D., Katsikis, V. N., & Stanimirovic, I. P. (2012). 
%            Full-rank representations of outer inverses based on the QR decomposition. 
%            Applied Mathematics and Computation, 218(20), 10321-10333.

  At=A';
  [Q,R,P]=qr(At);
  r=sum(any(abs(R)>1e-005,2));
  Q1=Q(:,1:r);
  R1=R(1:r,:);
  G=R1*P';
  Y=linsolve(G*A*Q1,G);
  X=Q1*Y;
end
