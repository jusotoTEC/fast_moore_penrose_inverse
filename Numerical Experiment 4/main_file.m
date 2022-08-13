function main_file()

    % Numerical Experiment 4

    % Reference: Soto-Quiros, P. (2022), A fast method to estimate the Moore-Penrose 
    %            inverse for well-determined numerical rank matrices based on the 
    %            Tikhonov regularization. (Submitted paper)

    dim_m=[20 100 200 300:300:6000];
    time_sec=[];
    time_par=[];


    for m=dim_m
        display(['Running experiment for m = ' num2str(m)])
        n=10000; r=floor(2*min([m n])/3);
        A = randn(m,r)*rand(r,n);
        ADist = distributed(A);

        b = sum(A,2);
        bDist = sum(ADist,2);

        xEx = ones(n,1);
        xDistEx = ones(n,1,'distributed');

        tic
        x = proposed_method(A,eps)*b;
        t1=toc;
        time_sec=[time_sec t1];

        tic
        xDist = mldivide(ADist,bDist);
        t2=toc;
        time_par=[time_par t2];
    end

    grid on
    hold on
    plot(dim_m,time_sec,'b')
    plot(dim_m,time_par,'r')
    xlabel('Dimension (m)')
    ylabel('Time (s)')
    legend('Alg. 1', 'Parallel Implementation')
end
