function main_file()

  % Numerical Experiment 1

  % Reference: Soto-Quiros, P. (2024), A fast method to estimate the Moore-Penrose
  %            inverse for well-determined numerical rank matrices based on the
  %            Tikhonov regularization. (Submitted paper)

  clc; clear; close all;
  warning ('off','all');
  tam = 5000:2500:20000;
  time=zeros(5,length(tam));
  error=zeros(4,length(tam));
  speed=zeros(4,length(tam));
  per_dif=zeros(4,length(tam));
  k=0;
  tol=eps;

  for m=tam

    k=k+1;

    %Create Matrix

    % Table 1 of paper
    r=round(m/4); A=randn(m,r)*randn(r,m/2);

    % Table 2 of paper
    %r=round(m/2); A=randn(m,r)*randn(r,m);

    % Table 3 of paper
    %A=randn(m,m/4);

    disp(['Running experiment for m = ', num2str(m)])

    %Time
    tic; A1=pinv(A); t1=toc; time(1,k)=t1;
    tic; A2=proposed_method(A,tol); t2=toc; time(2,k)=t2;
    tic; A3=geninv(A); t3=toc; time(3,k)=t3;
    tic; A4=qrginv(A); t4=toc; time(4,k)=t4;
    tic; A5=imqrginv(A); t5=toc; time(5,k)=t5;
    tic; A6=ats2(A); t6=toc; time(6,k)=t6;

    %Error
    error(1,k)=norm(A1-A2,'fro')^2;
    error(2,k)=norm(A1-A3,'fro')^2;
    error(3,k)=norm(A1-A4,'fro')^2;
    error(4,k)=norm(A1-A5,'fro')^2;
    error(5,k)=norm(A1-A6,'fro')^2;

    %Speed
    speed(1,k)=t1/t2;
    speed(2,k)=t1/t3;
    speed(3,k)=t1/t4;
    speed(4,k)=t1/t5;
    speed(5,k)=t1/t6;

    %Percentage Difference
    per_dif(1,k)=100*(1-t2/t1);
    per_dif(2,k)=100*(1-t3/t1);
    per_dif(3,k)=100*(1-t4/t1);
    per_dif(4,k)=100*(1-t5/t1);
    per_dif(5,k)=100*(1-t6/t1);
  end

  dim_tam=length(tam);
  dimension=cell(1,dim_tam);

  for k=1:dim_tam
      dimension{k}=num2str(tam(k));
  end

  %Table I =  Time Methods

  fprintf('Table 1: Time Methods\n')
  table_Results=table((time(1,:))',(time(2,:))',(time(3,:))',(time(4,:))', (time(5,:))',(time(6,:))', 'RowNames', dimension);
  table_Results.Properties.VariableNames={'Time_pinv', 'Time_pm', 'Time_geninv', 'Time_qrginv', 'Time_imqrginv','Time_ats2'};
  disp(table_Results)

  %Table II =  Error Methods

  fprintf('Table 2: Error Methods\n')
  table_Results=table((error(1,:))',(error(2,:))',(error(3,:))',(error(4,:))', (error(5,:))', 'RowNames', dimension);
  table_Results.Properties.VariableNames={'Error_pm', 'Error_geninv', 'Error_qrginv', 'Error_imqrginv','Error_ats2'};
  disp(table_Results)

  %Table III =  Speed Methods

  fprintf('Table 3: Speed Methods\n')
  table_Results=table((speed(1,:))',(speed(2,:))',(speed(3,:))',(speed(4,:))', (speed(5,:))', 'RowNames', dimension);
  table_Results.Properties.VariableNames={'Speed_pm', 'Speed_geninv', 'Speed_qrginv', 'Speed_imqrginv','Speed_ats2'};
  disp(table_Results)

  %Table IV =  Percent Difference Methods

  fprintf('Table 4: Percent Difference Methods\n')
  table_Results=table((per_dif(1,:))',(per_dif(2,:))',(per_dif(3,:))',(per_dif(4,:))', (per_dif(5,:))', 'RowNames', dimension);
  table_Results.Properties.VariableNames={'Per_Dif_pm', 'Per_Dif_geninv', 'Per_Dif_qrginv', 'Per_Dif_imqrginv','Per_Dif_ats2'};
  disp(table_Results)

end
