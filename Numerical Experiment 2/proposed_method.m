function X=proposed_method(A,s,tol)
% Method to compute Pseudoinverse of A
% Sintaxis:  X=proposed_method(A,svd,tol)
%              Input:  Matrix X of size m x n
%                      Positive singular values s of size r             
%                      Constante tol > 0
%              Output: Pseudoinverse X of size n x m
% Reference: Soto-Quiros, P. (2024), A fast method to estimate the Moore-Penrose 
% inverse for well-determined numerical rank matrices based on the Tikhonov 
% regularization. (Submitted paper)

    [m,n]=size(A);    
    At=A';
    r=length(s);
    if m>=n
      T=At*A;
      if r==n
          X=linsolve(T,At);
      else
          alpha=(0.5)*sqrt(tol/sum(1./s.^6));
          X=linsolve(T+alpha*eye(n),At);
      end
    else      
      T=A*At;
      if r==m
          X=At/T;
      else
          alpha=(0.5)*sqrt(tol/sum(1./s.^6));
          X=At/(T+alpha*eye(m));          
      end        
    end
end
