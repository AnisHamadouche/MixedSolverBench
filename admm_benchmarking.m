function admm_benchmarking
    clear;
    clc;
    close all;
    addpath './ADMM'
    
    %% Problem data
    
    % s = RandStream.create('mt19937ar','seed',0);
    % RandStream.setDefaultStream(s);
    scenario = 'R';
    if scenario == 'S'
        % synthetic data
        m = 100;       % number of examples
        n = 500;      % number of features
        
        x0 = randn(n,1);
        A = randn(m,n);
        A = A*spdiags(1./sqrt(sum(A.^2))',0,n,n); % normalize columns
        v = sqrt(0.001)*randn(m,1);
        b = A*x0 + v;
    else
        % real data
        m = 1000;
        noise_mean = 0;
        noise_var = 0.000001; % 10^{-5}
        inputImage = im2double(imread('./data/mi1helicopter.jpg'));
        inputImage = inputImage(:,:,1);
        % reduce resolution to fit in memory
        N = 1.5;
        [rows, columns, numColorChannels] = size(inputImage);
        numOutputRows = round(rows/N);
        numOutputColumns = round(columns/N);
        x0 = imresize(inputImage, [numOutputRows, numOutputColumns]);
        [x0m,x0n]=size(x0);
        x0 = x0(:);
        A = randn(m,size(x0,1));
        [m,n]=size(A);
        %g = imfilter(f,h,'conv','circular'); % blur
        b = imnoise(A*x0,'gaussian',noise_mean,noise_var); % ading noise
        imshow(reshape(x0,x0m,x0n));
    end
    
    fprintf('solving instance with %d examples, %d variables\n', m, n);
    fprintf('nnz(x0) = %d; signal-to-noise ratio: %.2f\n', nnz(x0), norm(A*x0)^2/noise_var);
    
    
    gamma_max = norm(A'*b,'inf');
    gamma = 0.1*gamma_max;
    
    % initialization
    
    x_init = rand(size(x0));
    z_init = x_init;
    u_init = x_init-z_init;
    
    % cached computations for all methods
    AtA = A'*A;
    Atb = A'*b;
    
    %% Global constants and defaults
    
    MAX_ITER = 500;
    ABSTOL   = 1e-4;
    RELTOL   = 1e-2;
    
    %% Extra parameters for ADMM
    lambda = 1;
    rho = 1/lambda;

    f = @(u) 0.5*sum_square(A*u-b);

    % Run 
    [y0,x64] = admm_entrypoint(A, b, Atb, lambda, gamma, rho, MAX_ITER, m, n, ABSTOL, RELTOL,'double');
    fprintf('Solved with double precision\n');
    [psnr_64,snr_64] =  psnr(x64,x0);
    ssimval_64 = ssim(x64,x0);  
    
    [y32,x32] = admm_entrypoint(A, b, Atb, lambda, gamma, rho, MAX_ITER, m, n, ABSTOL, RELTOL,'single');
    fprintf('Solved with single precision\n');
    [psnr_32,snr_32] =  psnr(x32,x0);
    ssimval_32 = ssim(x32,x0);  
    
    [y16,x16] = admm_entrypoint(x_init,z_init,u_init,A, b, Atb, lambda, gamma, MAX_ITER, m, n, ABSTOL, RELTOL,'fixed16');
    fprintf('Solved with half precision\n');
    [psnr_16,snr_16] =  psnr(x16,x0);
    ssimval_16 = ssim(x16,x0);
    
    [y12,x12] = admm_entrypoint(x_init,z_init,u_init,A, b, Atb, lambda, gamma, MAX_ITER, m, n, ABSTOL, RELTOL,'fixed12');
    fprintf('Solved with 12 bits\n');
    [psnr_12,snr_12] =  psnr(x12,x0);
    ssimval_12 = ssim(x12,x0);    
    
    [y8,x8] = admm_entrypoint(x_init,z_init,u_init,A, b, Atb, lambda, gamma, MAX_ITER, m, n, ABSTOL, RELTOL,'fixed8');
    fprintf('Solved with 8 bits\n');
    [psnr_8,snr_8] =  psnr(x8,x0);
    ssimval_8 = ssim(x8,x0);
    
    [y2,x2] = admm_entrypoint(x_init,z_init,u_init,A, b, Atb, lambda, gamma, MAX_ITER, m, n, ABSTOL, RELTOL,'fixed2');
    fprintf('Solved with 2 bits\n');
    [psnr_2,snr_2] =  psnr(x2,x0);
    ssimval_2 = ssim(x2,x0);
    
    figure;
    subplot(5,1,1);semilogy(y0,'k'); 
    legend('Ground truth objective') 
    title('Ground truth (Double precision)') 
    
    subplot(5,2,3);semilogy(y0,'k'); 
    title('Single precision floating-point output') 
    subplot(5,2,4);
    try
      semilogy(abs(y0(1:numel(y32))-double(y32)),'r');
      title('Residual error');
    catch
      semilogy(abs(y0-double(y32(1:numel(y0)))),'r');
      title('Residual error');
    end
    
    subplot(5,2,5);semilogy(y16,'k'); 
    title('16-bit fixed-point output') 
    subplot(5,2,6);
    try
      semilogy(abs(y0(1:numel(y16))-double(y16)),'r');
    catch
      semilogy(abs(y0-double(y16(1:numel(y0)))),'r');
    end
    
    subplot(5,2,7);semilogy(y8,'k'); 
    title('8-bit fixed-point output') 
    subplot(5,2,8);
    try
      semilogy(abs(y0(1:numel(y8))-double(y8)),'r');
    catch
      semilogy(abs(y0-double(y8(1:numel(y0)))),'r');
    end
    
    subplot(5,2,9);semilogy(y2,'k'); 
    title('2-bit fixed-point output') 
    subplot(5,2,10);
    try
      semilogy(abs(y0(1:numel(y2))-double(y2)),'r');
    catch
      semilogy(abs(y0-double(y2(1:numel(y0)))),'r');
    end
    
    xlabel('Iterations, k') 
    
    figure;
    %imshow([reshape(x0,x0m,x0n),reshape(x64,x0m,x0n), ...
    %    reshape(x32,x0m,x0n),reshape(x16,x0m,x0n),reshape(x12,x0m,x0n),reshape(x8,x0m,x0n)],'InitialMagnification', 800);
    h64 = subplot(2,3,1), imshow(reshape(x64,x0m,x0n), 'Parent', h64), title(h64, ['BW = 64 bits (',num2str(round(100*64/64)),'%); SSIM = ', num2str(ssimval_64),'; PSNR = ', num2str(psnr_64)]); hold  on
    h32 = subplot(2,3,2), imshow(reshape(x32,x0m,x0n), 'Parent', h32), title(h32, ['BW = 32 bits (',num2str(round(100*32/64)),'%); SSIM = ', num2str(ssimval_32),'; PSNR = ', num2str(psnr_32)]); hold  on
    h16 = subplot(2,3,3), imshow(reshape(x16,x0m,x0n), 'Parent', h16), title(h16, ['BW = 16 bits (',num2str(round(100*16/64)),'%); SSIM = ', num2str(ssimval_16),'; PSNR = ', num2str(psnr_16)]); hold  on
    h12 = subplot(2,3,4), imshow(reshape(x12,x0m,x0n), 'Parent', h12), title(h12, ['BW = 12 bits (',num2str(round(100*12/64)),'%); SSIM = ', num2str(ssimval_12),'; PSNR = ', num2str(psnr_12)]); hold  on
    h8 = subplot(2,3,5), imshow(reshape(x8,x0m,x0n), 'Parent', h8), title(h8, ['BW = 8 bits (',num2str(round(100*8/64)),'%); SSIM = ', num2str(ssimval_8),'; PSNR = ', num2str(psnr_8)]); hold  on
    h2 = subplot(2,3,6), imshow(reshape(x2,x0m,x0n), 'Parent', h2), title(h2, ['BW = 2 bits (',num2str(round(100*2/64)),'%); SSIM = ', num2str(ssimval_2),'; PSNR = ', num2str(psnr_2)]); hold  off
    
    
    admm_solv_mex
end 