function X = geninv(A)
% Method to compute Pseudoinverse of A
% Sintaxis:  X = geninv(A)
%                Input:  Matrix X of size m x n (m>=n)
%                Output: Pseudoinverse X of size n x m
% Reference: Courrieu, P. (2005). Fast Computation of Moore-Penrose Inverse 
%            Matrices. Neural Information Processing-Letters and Reviews, 8(2).
    n=size(A,2);
    G=A'*A;
    dG=diag(G);
    tol= min(dG(dG>0))*1e-9;
    L=zeros(size(G));
    r=0;
    for k=1:n
        r=r+1;
        L(k:n,r)=G(k:n,k)-L(k:n,1:(r-1))*L(k,1:(r-1))';
        if L(k,r)>tol
            L(k,r)=sqrt(L(k,r));
            if k<n
                L((k+1):n,r)=L((k+1):n,r)/L(k,r);
            end
        else
            r=r-1;
        end
    end
    L1=L(:,1:r);
    M=linsolve(L1'*L1,eye(r));
    X=L1*(M*M)*L1'*A';
end


