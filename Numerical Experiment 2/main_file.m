function main_file()

  % Numerical Experiment 2
  
  % Reference: Soto-Quiros, P. (2024), A fast method to estimate the Moore-Penrose 
  %            inverse for well-determined numerical rank matrices based on the 
  %            Tikhonov regularization. (Submitted paper)

  clc; clear; close all;
  warning ('off','all');
  tam = 1000:1000:10000;
  time=zeros(5,length(tam));
  error=zeros(4,length(tam));
  speed=zeros(4,length(tam));
  per_dif=zeros(4,length(tam));
  k=0;
  tol=eps;
  for m=tam
      
    k=k+1;
    
    %Create Matrix and SVD
    
    [A,s]=structured_matrix(m);
    
    disp(['Running experiment for m = ', num2str(m)])

    %Time
    tic; A1=pinv(A); t1=toc; time(1,k)=t1;
    tic; A2=proposed_method(A,s,tol); t2=toc; time(2,k)=t2;
    tic; A3=qrginv(A); t3=toc; time(3,k)=t3;
    tic; A4=imqrginv(A); t4=toc; time(4,k)=t4;
    tic; A5=ats2(A); t5=toc; time(5,k)=t5;
    
    %Error
    error(1,k)=norm(A1-A2,'fro')^2;
    if error(1,k)>eps
        disp('VOLVER A CORRER')
        break
    end
    error(2,k)=norm(A1-A3,'fro')^2;
    error(3,k)=norm(A1-A4,'fro')^2;
    error(4,k)=norm(A1-A5,'fro')^2;
    
    %Speed
    speed(1,k)=t1/t2;
    speed(2,k)=t1/t3;
    speed(3,k)=t1/t4;
    speed(4,k)=t1/t5;
    
    %Percentage Difference
    per_dif(1,k)=100*(1-t2/t1);
    per_dif(2,k)=100*(1-t3/t1);
    per_dif(3,k)=100*(1-t4/t1);
    per_dif(4,k)=100*(1-t5/t1);
  end
  
  dim_tam=length(tam);
  dimension=cell(1,dim_tam);
  
  for k=1:dim_tam
      dimension{k}=num2str(tam(k));
  end
    
  %Table I =  Time Methods
  
  fprintf('Table 1: Time Methods\n')
  table_Results=table((time(1,:))',(time(2,:))',(time(3,:))',(time(4,:))', (time(5,:))', 'RowNames', dimension);
  table_Results.Properties.VariableNames={'Time_pinv', 'Time_pm', 'Time_qrginv', 'Time_imqrginv','Time_ats2'};
  disp(table_Results)
  
  %Table II =  Error Methods
  
  fprintf('Table 2: Error Methods\n')
  table_Results=table((error(1,:))',(error(2,:))',(error(3,:))',(error(4,:))',  'RowNames', dimension);
  table_Results.Properties.VariableNames={'Error_pm', 'Error_qrginv', 'Error_imqrginv','Error_ats2'};
  disp(table_Results)
  
  %Table III =  Speed Methods
  
  fprintf('Table 3: Speed Methods\n')
  table_Results=table((speed(1,:))',(speed(2,:))',(speed(3,:))',(speed(4,:))', 'RowNames', dimension);
  table_Results.Properties.VariableNames={'Speed_pm', 'Speed_qrginv', 'Speed_imqrginv','Speed_ats2'};
  disp(table_Results)  
  
  %Table IV =  Percent Difference Methods
  
  fprintf('Table 4: Percent Difference Methods\n')
  table_Results=table((per_dif(1,:))',(per_dif(2,:))',(per_dif(3,:))',(per_dif(4,:))', 'RowNames', dimension);
  table_Results.Properties.VariableNames={'Per_Dif_pm', 'Per_Dif_qrginv', 'Per_Dif_imqrginv','Per_Dif_ats2'};
  disp(table_Results)  
  
  
  %Diagrams
  
  %Figure 1 - Time
  figure
  hold on  
  for k=1:5
    plot(tam,time(k,:))    
  end
  legend('pinv','Alg. 1','qrginv','imqrginv','ats2')
  xlabel('Dimension (m)')
  ylabel('Time (s)')
  grid on
  pbaspect([2 1 1])
  
  %Figure 2 - Error
  figure
  hold on  
  for k=1:4
    plot(tam,error(k,:))    
  end
  legend('Alg. 1','qrginv','imqrginv','ats2')
  xlabel('Dimension (m)')
  ylabel('Error')
  grid on
  pbaspect([2 1 1])
  
  %Figure 3 - Speed
  figure
  hold on  
  for k=1:4
    plot(tam,speed(k,:))    
  end
  legend('Alg. 1','qrginv','imqrginv','ats2')
  xlabel('Dimension (m)')
  ylabel('Speed')
  grid on
  pbaspect([2 1 1])
  
  %Figure 4 - Error
  figure
  hold on  
  for k=1:4
    plot(tam,per_dif(k,:))    
  end
  legend('Alg. 1','qrginv','imqrginv','ats2')
  xlabel('Dimension (m)')
  ylabel('Percent Difference')
  grid on
  pbaspect([2 1 1])
    
end

function [A,s]=structured_matrix(m)    
    %Create structured matrix given in equation (23) of the paper
    d=100*rand(1);
    d0=[1 (d+1)*ones(1,m-2) d];
    dp=sqrt(d)*ones(m-1,1);
    A=diag(dp,-1)+diag(d0,0)+diag(dp,1);
    v=1:m-1;
    s=sort(abs(d+1+2*sqrt(d)*cos(pi*v/m)),'descend');
end

