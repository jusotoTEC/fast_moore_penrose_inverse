function X = imqrginv(A)
% Method to compute Pseudoinverse of A
% Sintaxis:  X = imqrginv(A)
%                Input:  Matrix X of size m x n (m>=n)
%                Output: Pseudoinverse X of size n x m
% Reference: Ataei, A. (2014). Improved Qrginv algorithm for computing Moore-Penrose 
%            inverse matrices. International Scholarly Research Notices, 2014.

  [Q,R,P]=qr(A);
  r=sum(any(abs(R)>1e-005,2));
  Q1=Q(:,1:r);
  R1=R(1:r,:);
  Rt=R1';
  S=linsolve(R1*Rt,Q1');
  X=P*Rt*S;
end
