function X = qrginv(A)
% Method to compute Pseudoinverse of A
% Sintaxis:  X = urging(A)
%                Input:  Matrix X of size m x n (m>=n)
%                Output: Pseudoinverse X of size n x m
% Reference: Katsikis, V. N., Pappas, D., & Petralias, A. (2011). An improved 
%            method for the computation of the Moore-Penrose inverse matrix. 
%            Applied Mathematics and Computation, 217(23), 9828-9834.

    [m,n] = size(A);
    [Q,R,P] = qr(A);
    r = sum(any(abs(R)>1e-005,2));
    R1 = R(1:r,:);
    C = R1*R1';
    G=linsolve(C,R1);
    R3 = [G' zeros(n,m-r)];
    X = P*R3*Q';
end
