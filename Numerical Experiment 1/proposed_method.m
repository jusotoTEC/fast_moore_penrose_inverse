function X=proposed_method(A,tol)
% Method to compute Pseudoinverse of A
% Sintaxis:  X=proposed_method(A,tol)
%              Input:  Matrix X of size m x n
%                      Constante tol > 0
%              Output: Pseudoinverse X of size n x m
% Reference: Soto-Quiros, P. (2024), A fast method to estimate the Moore-Penrose 
% inverse for well-determined numerical rank matrices based on the Tikhonov 
% regularization. 

    [m,n]=size(A);    
    At=A';
    if m>=n
      T=At*A;
      num_cond=rcond(T);
      if num_cond>sqrt(eps)
          X=linsolve(T,At);
      else
          s1=eig(T);
          s2=sort(s1,'descend');
          f=n*eps*s2(1);
          s=sqrt(s1(s1>f));         
          alpha=(0.5)*sqrt(tol/sum(1./s.^6));
          X=linsolve(T+alpha*eye(n),At);
      end
    else      
      T=A*At;
      num_cond=rcond(T);
      if num_cond>sqrt(eps)
          X=At/T;
      else
          s1=eig(T);
          s2=sort(s1,'descend');
          f=m*eps*s2(1);
          s=sqrt(s2(s2>f));         
          alpha=(0.5)*sqrt(tol/sum(1./s.^6));
          X=At/(T+alpha*eye(m));          
      end        
    end
end
