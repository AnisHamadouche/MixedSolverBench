function admm_benchmark
    clear;
    clc;
    close all;
    addpath './ADMM'

    %% configuraion
    report = false; % set to 'true' for instrumentation report
    scenario = 'S'; % 'R': real, 'S': synthetic data
    
    %% Problem data
    
    % s = RandStream.create('mt19937ar','seed',0);
    % RandStream.setDefaultStream(s);
    
    if scenario == 'S'
        % synthetic data
        n = 500;      % number of features
        noise_var=0.001;
        
        % sparse vector with 0.1*n nonzero entries
        x0 = sprand(n,1,0.1);
        s0 = 0.1*n;
        %m > 2*s0*log(n/s0) + (7/5)*s0 + 1
        m = round(2*s0*log(n/s0) + (7/5)*s0 + 1)+10; % number of examples
        
        H = randn(m,n);
        H = H*spdiags(1./sqrt(sum(H.^2))',0,n,n); % normalize columns
        v = sqrt(noise_var)*randn(m,1);
        b = H*x0 + v;
    else
        % real data
        noise_mean = 0;
        noise_var = 0.000001; % 10^{-5}
        inputImage = im2double(imread('./data/mi1helicopter.jpg'));
        inputImage = inputImage(:,:,1);
        % reduce resolution to fit in memory
        N = 3;
        [rows, columns, numColorChannels] = size(inputImage);
        numOutputRows = round(rows/N);
        numOutputColumns = round(columns/N);
        x0 = imresize(inputImage, [numOutputRows, numOutputColumns]);
        [x0m,x0n]=size(x0);
        x0 = x0(:);
        s0 = 0.1*size(x0,1);
        m = round(2*s0*log(size(x0,1)/s0) + (7/5)*s0 + 1)+10; % number of examples
        H = round(rand(m,size(x0,1)));
        %A = eye(size(x0,1));
        [m,n]=size(H);
        %g = imfilter(f,h,'conv','circular'); % blur
        %b = imnoise(A*x0,'gaussian',noise_mean,noise_var); % ading noise
        v = sqrt(0.001)*randn(m,1);
        b = H*x0 + v;
        imshow(reshape(x0,x0m,x0n));
    end
    
    fprintf('solving instance with %d examples, %d variables\n', m, n);
    fprintf('nnz(x0) = %d; signal-to-noise ratio: %.2f\n', nnz(x0), norm(H*x0)^2/noise_var);
    
    
    gamma_max = norm(H'*b,'inf');
    gamma = 0.1*gamma_max;
    
    % initialization
    
    x_init = rand(size(x0));
    z_init = x_init;
    u_init = x_init-z_init;
    
    % cached computations for all methods
    AtA = H'*H;
    Atb = H'*b;
    
    %% Global constants and defaults
    
    MAX_ITER = 1000;
    ABSTOL   = 1e-4;
    RELTOL   = 1e-2;
    
    %% Extra parameters for ADMM
    lambda = 1;
    rho = 1/lambda;

    % Run 

    [y64,x64,z64,u64,y64_ref,x64_ref,z64_ref,u64_ref] = admml0_entrypoint(H, b, Atb, lambda, gamma, rho, MAX_ITER, m, n, ABSTOL, RELTOL,'double');
    save('./data/admm/admm_xk_fp64.mat', 'x64');
    save('./data/admm/admm_zk_fp64.mat', 'z64');
    save('./data/admm/admm_uk_fp64.mat', 'u64');
    save('./data/admm/admm_fk_fp64.mat', 'y64');
    save('./data/admm/admm_xk_ref_fp64.mat', 'x64_ref');
    save('./data/admm/admm_zk_ref_fp64.mat', 'z64_ref');
    save('./data/admm/admm_uk_ref_fp64.mat', 'u64_ref');
    save('./data/admm/admm_fk_ref_fp64.mat', 'y64_ref');  
    fprintf('Solved with double precision\n');

    [y32,x32,z32,u32,y32_ref,x32_ref,z32_ref,u32_ref] = admml0_entrypoint(H, b, Atb, lambda, gamma, rho, MAX_ITER, m, n, ABSTOL, RELTOL,'single');
    save('./data/admm/admm_xk_fp32.mat', 'x32');
    save('./data/admm/admm_zk_fp32.mat', 'z32');
    save('./data/admm/admm_uk_fp32.mat', 'u32');
    save('./data/admm/admm_fk_fp32.mat', 'y32');
    save('./data/admm/admm_xk_ref_fp32.mat', 'x32_ref');
    save('./data/admm/admm_zk_ref_fp32.mat', 'z32_ref');
    save('./data/admm/admm_uk_ref_fp32.mat', 'u32_ref');
    save('./data/admm/admm_fk_ref_fp32.mat', 'y32_ref');   
    fprintf('Solved with single precision\n');

    [y16,x16,z16,u16,y16_ref,x16_ref,z16_ref,u16_ref] = admml0_entrypoint(H, b, Atb, lambda, gamma, rho, MAX_ITER, m, n, ABSTOL, RELTOL,'fixed16');
    save('./data/admm/admm_xk_fx16.mat', 'x16');
    save('./data/admm/admm_zk_fx16.mat', 'z16');
    save('./data/admm/admm_uk_fx16.mat', 'u16');
    save('./data/admm/admm_fk_fx16.mat', 'y16');
    save('./data/admm/admm_xk_ref_fx16.mat', 'x16_ref');
    save('./data/admm/admm_zk_ref_fx16.mat', 'z16_ref');
    save('./data/admm/admm_uk_ref_fx16.mat', 'u16_ref');
    save('./data/admm/admm_fk_ref_fx16.mat', 'y16_ref');
    fprintf('Solved with half precision\n');

    [y12,x12,z12,u12,y12_ref,x12_ref,z12_ref,u12_ref] = admml0_entrypoint(H, b, Atb, lambda, gamma, rho, MAX_ITER, m, n, ABSTOL, RELTOL,'fixed12');
    save('./data/admm/admm_xk_fx12.mat', 'x12');
    save('./data/admm/admm_zk_fx12.mat', 'z12');
    save('./data/admm/admm_uk_fx12.mat', 'u12');
    save('./data/admm/admm_fk_fx12.mat', 'y12'); 
    save('./data/admm/admm_xk_ref_fx12.mat', 'x12_ref');
    save('./data/admm/admm_zk_ref_fx12.mat', 'z12_ref');
    save('./data/admm/admm_uk_ref_fx12.mat', 'u12_ref');
    save('./data/admm/admm_fk_ref_fx12.mat', 'y12_ref'); 
    fprintf('Solved with 12 bits\n');    
    
    [y8,x8,z8,u8,y8_ref,x8_ref,z8_ref,u8_ref] = admml0_entrypoint(H, b, Atb, lambda, gamma, rho, MAX_ITER, m, n, ABSTOL, RELTOL,'fixed8');
    save('./data/admm/admm_xk_fx8.mat', 'x8');
    save('./data/admm/admm_zk_fx8.mat', 'z8');
    save('./data/admm/admm_uk_fx8.mat', 'u8');
    save('./data/admm/admm_fk_fx8.mat', 'y8');    
    save('./data/admm/admm_xk_ref_fx8.mat', 'x8_ref');
    save('./data/admm/admm_zk_ref_fx8.mat', 'z8_ref');
    save('./data/admm/admm_uk_ref_fx8.mat', 'u8_ref');
    save('./data/admm/admm_fk_ref_fx8.mat', 'y8_ref'); 
    fprintf('Solved with 8 bits\n');
    
    [y2,x2,z2,u2,y2_ref,x2_ref,z2_ref,u2_ref] = admml0_entrypoint(H, b, Atb, lambda, gamma, rho, MAX_ITER, m, n, ABSTOL, RELTOL,'fixed2');
    save('./data/admm/admm_xk_fx2.mat', 'x2');
    save('./data/admm/admm_zk_fx2.mat', 'z2');
    save('./data/admm/admm_uk_fx2.mat', 'u2');
    save('./data/admm/admm_fk_fx2.mat', 'y2');
    save('./data/admm/admm_xk_ref_fx2.mat', 'x2_ref');
    save('./data/admm/admm_zk_ref_fx2.mat', 'z2_ref');
    save('./data/admm/admm_uk_ref_fx2.mat', 'u2_ref');
    save('./data/admm/admm_fk_ref_fx2.mat', 'y2_ref');
    fprintf('Solved with 2 bits\n');

    if (strcmp(scenario,'R'))
        [psnr_64,snr_64] =  psnr(x64(:,end),x0)
        ssimval_64 = ssim(x64(:,end),x0)
        imshow(reshape(x64(:,end),x0m,x0n));
        [psnr_32,snr_32] =  psnr(x32(:,end),x0);
        ssimval_32 = ssim(x32(:,end),x0); 
        [psnr_16,snr_16] =  psnr(x16(:,end),x0);
        ssimval_16 = ssim(x16(:,end),x0);
        [psnr_12,snr_12] =  psnr(x12(:,end),x0);
        ssimval_12 = ssim(x12(:,end),x0);
        [psnr_8,snr_8] =  psnr(x8(:,end),x0);
        ssimval_8 = ssim(x8(:,end),x0);
        [psnr_2,snr_2] =  psnr(x2(:,end),x0);
        ssimval_2 = ssim(x2(:,end),x0);
        figure;
        %imshow([reshape(x0,x0m,x0n),reshape(x64,x0m,x0n), ...
        %    reshape(x32,x0m,x0n),reshape(x16,x0m,x0n),reshape(x12,x0m,x0n),reshape(x8,x0m,x0n)],'InitialMagnification', 800);
        h64 = subplot(2,3,1), imshow(reshape(x64,x0m,x0n), 'Parent', h64), title(h64, ['BW = 64 bits (',num2str(round(100*64/64)),'%); SSIM = ', num2str(ssimval_64),'; PSNR = ', num2str(psnr_64)]); hold  on
        h32 = subplot(2,3,2), imshow(reshape(x32,x0m,x0n), 'Parent', h32), title(h32, ['BW = 32 bits (',num2str(round(100*32/64)),'%); SSIM = ', num2str(ssimval_32),'; PSNR = ', num2str(psnr_32)]); hold  on
        h16 = subplot(2,3,3), imshow(reshape(x16,x0m,x0n), 'Parent', h16), title(h16, ['BW = 16 bits (',num2str(round(100*16/64)),'%); SSIM = ', num2str(ssimval_16),'; PSNR = ', num2str(psnr_16)]); hold  on
        h12 = subplot(2,3,4), imshow(reshape(x12,x0m,x0n), 'Parent', h12), title(h12, ['BW = 12 bits (',num2str(round(100*12/64)),'%); SSIM = ', num2str(ssimval_12),'; PSNR = ', num2str(psnr_12)]); hold  on
        h8 = subplot(2,3,5), imshow(reshape(x8,x0m,x0n), 'Parent', h8), title(h8, ['BW = 8 bits (',num2str(round(100*8/64)),'%); SSIM = ', num2str(ssimval_8),'; PSNR = ', num2str(psnr_8)]); hold  on
        h2 = subplot(2,3,6), imshow(reshape(x2,x0m,x0n), 'Parent', h2), title(h2, ['BW = 2 bits (',num2str(round(100*2/64)),'%); SSIM = ', num2str(ssimval_2),'; PSNR = ', num2str(psnr_2)]); hold  off
    end
    
    figure;
    subplot(5,1,1);semilogy(y64,'k'); 
    legend('Ground truth objective') 
    title('Ground truth (Double precision)') 
    
    subplot(5,2,3);semilogy(y32,'k'); 
    title('Single precision floating-point output') 
    subplot(5,2,4);
    try
      semilogy(abs(y64(1:numel(y32))-double(y32)),'r');
      title('Residual error');
    catch
      semilogy(abs(y64-double(y32(1:numel(y64)))),'r');
      title('Residual error');
    end
    
    subplot(5,2,5);semilogy(y16,'k'); 
    title('16-bit fixed-point output') 
    subplot(5,2,6);
    try
      semilogy(abs(y64(1:numel(y16))-double(y16)),'r');
    catch
      semilogy(abs(y64-double(y16(1:numel(y64)))),'r');
    end
    
    subplot(5,2,7);semilogy(y8,'k'); 
    title('8-bit fixed-point output') 
    subplot(5,2,8);
    try
      semilogy(abs(y64(1:numel(y8))-double(y8)),'r');
    catch
      semilogy(abs(y64-double(y8(1:numel(y64)))),'r');
    end
    
    subplot(5,2,9);semilogy(y2,'k'); 
    title('2-bit fixed-point output') 
    subplot(5,2,10);
    try
      semilogy(abs(y64(1:numel(y2))-double(y2)),'r');
    catch
      semilogy(abs(y64-double(y2(1:numel(y64)))),'r');
    end
    
    xlabel('Iterations, k') 
    
    % reporting
    if(report)
        admm_solv_mex
    end

    para.gamma = gamma;
    para.rho = rho;
    verify('fp64',para)
    verify('fp32',para)
    verify('fp16',para)
    verify('fp12',para)
    verify('fp8',para)
end 
