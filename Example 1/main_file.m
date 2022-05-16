% Example 1

% Reference: Soto-Quiros, P. (2022), A fast method to estimate the Moore-Penrose 
%            inverse for well-determined numerical rank matrices based on the 
%            Tikhonov regularization. (Submitted paper)

clc; clear; 

%Unknown vector bt
bt=[0.53 0.97 1.06 0.40 1.20]';

% Information available
A=[8 10 19 16;
  31 26 12 28;
  16 20 38 32;
   7  8 13 12;
  21 24 39 36];
er=[0.02168 0.08861 0.11303 0.11678 0.10061]';
b=bt+er;
tol=10^-4;

% Proposed Solution
B=sym(A'*A);
s_v=sort(sqrt(eig(B)),'descend');
s=s_v(1:2);
alpha1=double((norm(er)/norm(b))*sqrt(tol/norm(er)^2- sum(1./s.^2))/sqrt(sum(1./s.^6)));
alpha2=(floor(alpha1*10000)/10000)/2;
x2=(A'*A+alpha2*eye(4))\(A'*b);

% Solution of free-noisy system with pinv command
x1=pinv(A)*bt;

% Error
e1=norm(x1-x2)^2;

% Numerical verification
cond=e1<tol
