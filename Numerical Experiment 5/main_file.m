function main_file()

  % Reference: Soto-Quiros, P. (2024), A fast method to estimate the Moore-Penrose 
  %            inverse for well-determined numerical rank matrices based on the 
  %            Tikhonov regularization. (Submitted paper)

    clc; clear;
    
    %Matrix 
    A=[-1 1 -1 -2; 0 0 0 4; 2 -2 2 0; 0 0 0 -2; 1 -1 1 0];

    % Rank of matrix A
    r=rank(A);

    % Positive Singular Values of Matrix A
    aux=svd(A); sExact=aux(1:r);


    
    %Compute singular values using QR method
    T=A.'*A;    
    Tk=T;
    Uk=eye(size(T,1));
    iterMax=100;
    for k=1:iterMax
      [Qk,Rk]=qr(Tk);
      Tk=Rk*Qk;
      s1=diag(Tk);
      Uk=Uk*Qk;      
      s2=sort(s1,'descend');
      sk=sqrt(s2(1:r));      %Estimation of singular value of A                    
      % Verify condition
      condS=sum((sExact.^6-sk.^6)./((sk.^6).*(sExact.^6)));
      if condS>0
          break
      end   
    end

    tol=1e-10;
    alpha=(0.5)*sqrt(tol/sum(1./sk.^6));
    
    %Approximation Moore-Penrose with approximation of singular values
    %using QR method
    Xp1=linsolve(T+alpha*eye(4),A.');   
    
    %Approximation Moore-Penrose using singular values of eig command
    Xp2=proposed_method(A,tol);
  
    erQR=norm(Xp1-pinv(A),'fro')^2;
    
    er2SVD=norm(Xp2-pinv(A),'fro')^2;
    
    disp(['Error using approximation of singular values: ', num2str(erQR)])
    disp(['Error using exact values of singular values: ', num2str(er2SVD)])
end
