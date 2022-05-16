function main_file()

    % Numerical Experiment 3

    % Reference: Soto-Quiros, P. (2022), A fast method to estimate the Moore-Penrose 
    %            inverse for well-determined numerical rank matrices based on the 
    %            Tikhonov regularization. (Submitted paper)

    clc; clear;  close all   
    numImg=1110; %Number of images in folders "image_free_noise" and "image_with_noise"
    
    % Load training images    
    X=zeros(128^2,numImg); % Images Free-Noisy
    Y=zeros(128^2,numImg); % Image with Noise

    for k=1:numImg
        X1=imread(['image_free_noise/coast (' num2str(k) ').jpg']);
        X2=im2double(X1);    
        X(:,k)=X2(:);
        Y1=imread(['image_with_noise/coast (' num2str(k) ').jpg']);
        Y2=im2double(Y1);    
        Y(:,k)=Y2(:);        
    end
    
    numS=randi([1 numImg], [1 6]);
    
    %Show random images (free-noise)    
    figure        
    for k=1:length(numS)
       subplot(1,length(numS),k)
       imshow(reshape(X(:,numS(k)),[128 128]))
    end
    
    %Show random images (with noise)    
    figure    
    for k=1:length(numS)
       subplot(1,length(numS),k)
       imshow(reshape(Y(:,numS(k)),[128 128]))
    end
    
    % Check Matrix T is rank-deficient    
    T=Y'*Y;
    [m,n]=size(T);
    disp(['Dimension of matrix T = ' num2str(m) ' x ' num2str(n)])
    r=rank(T);
    disp(['Matrix T is rank-deficient because rank(T) = ' num2str(r)])

    %Compute Filter F1 (using proposed method) and F2 (using pinv)

    tic; Yp1=proposed_method(Y,eps); t1=toc;
    F1=X*Yp1; 
    disp(['Execution time to compute Moore-Penrose of Y using proposed_method = ' num2str(t1,8) ' seconds'])

    tic; Yp2=pinv(Y); t2=toc;
    disp(['Execution time to compute Moore-Penrose of Y using pinv = ' num2str(t2,8) ' seconds'])
    F2=X*Yp2;

    %Performe analysis of proposed_method vrs command pinv
    speedup=t2/t1;
    per_dif=100*(t2-t1)/t2;
    disp(['Speedup to compute Moore-Penrose of Y using proposed_method = ' num2str(speedup,8) ' seconds'])
    disp(['proposed_method is ' num2str(per_dif,8) '% faster than command pinv'])
    
    %Error
    error=norm(Yp1-Yp2)^2;
    disp(['Error to compute Moore-Penrose of Y using proposed_method = ' num2str(error)])

    %Clean noisy images - There are four test images   
    
    A1=imread('test_images\test_image (1).jpg');
    %A1=imread('test_images\test_image (2).jpg');
    %A1=imread('test_images\test_image (3).jpg');
    %A1=imread('test_images\test_image (4).jpg');
    
    Xt=im2double(A1);
    figure
    imshow(Xt)
    title('Source Image')
    Yt=Xt+0.1*randn(size(Xt));
    figure
    imshow(im2uint8(Yt))
    title('Noisy Image')
    yt_v=Yt(:);

    %proposed_method
    xt_v_pm=F1*yt_v;
    Xt_est_pm=im2uint8(reshape(xt_v_pm,size(Xt)));
    figure
    imshow(Xt_est_pm)
    title('Estimate Image with proposed\_method')
    error_estimation_pm=norm(im2double(Xt)-im2double(Xt_est_pm),'fro')/norm(im2double(Xt),'fro');
    disp(['Error estimation of proposed_method = ' num2str(error_estimation_pm,8)])

    %Comand pinv
    xt_v_pinv=F2*yt_v;
    Xt_est_pinv=im2uint8(reshape(xt_v_pinv,size(Xt)));
    figure
    imshow(Xt_est_pinv)
    title('Estimate Image with pinv')
    error_estimation_pinv=norm(im2double(Xt)-im2double(Xt_est_pinv),'fro')/norm(im2double(Xt),'fro');
    disp(['Error estimation of pinv = ' num2str(error_estimation_pinv,8)])
                           
end