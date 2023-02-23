%Use pg_benchmarking with custom problem data (modify the data) or randomly 
%generated problem data (default) to compare the LASSO function values
%under different machine representation. The default function implements PG
%under 'double precision', 'single precision', '12 bits fixed-point' and 
%'16 bits fixed-point' representations. Edit 'mytypes.m' to add custom data
%types. type 'help mytypes.m' for furhter information.

function pg_benchmarking
clear;
clc;
close all;
addpath './PG'
addpath './utils'

%% Problem data

% s = RandStream.create('mt19937ar','seed',0);
% RandStream.setDefaultStream(s);
scenario = 'R';
if (scenario == 'S')
    % synthetic data
    n = 500;      % number of features
    noise_var=0.001;
    
    % sparse vector with 0.1*n nonzero entries
    x0 = sprand(n,1,0.1);
    s0 = 0.1*n;
    %m > 2*s0*log(n/s0) + (7/5)*s0 + 1
    m = round(2*s0*log(n/s0) + (7/5)*s0 + 1)+10; % number of examples
    A = randn(m,n);
    A = A*spdiags(1./sqrt(sum(A.^2))',0,n,n); % normalize columns
    v = sqrt(noise_var)*randn(m,1);
    b = A*x0 + v;
else
    % real data
    noise_mean = 0;
    noise_var = 0.0001; % 10^{-5}
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
    n = size(x0,1);
    s0 = 0.0001*n;
    m = round(2*s0*log(n/s0) + (7/5)*s0 + 1)+10; % number of examples
    A = randn(m,n);
    [m,n]=size(A);
    %g = imfilter(f,h,'conv','circular'); % blur
    b = imnoise(A*x0,'gaussian',noise_mean,noise_var); % ading noise
    %imshow(reshape(x0,x0m,x0n));
end

fprintf('solving instance with %d examples, %d variables\n', m, n);
fprintf('nnz(x0) = %d; signal-to-noise ratio: %.2f\n', nnz(x0), norm(A*x0)^2/noise_var);

gamma_max = norm(A'*b,'inf');
gamma = 0.1*gamma_max;

% cached computations for all methods
AtA = A'*A;
Atb = A'*b;

%% Global constants and defaults

MAX_ITER = 100;
ABSTOL   = 1e-4;
RELTOL   = 1e-2;

%% Extra parameters for ADMM
lambda = 1;
rho = 1/lambda;

%   % Test inputs  
%   m = 200;       % number of examples
%   n = 100;      % number of features
%   %x0 = sprandn(n,1,0.05);
%   x0 = randn(n,1);
%   A = randn(m,n);
%   %A = A*spdiags(1./sqrt(sum(A.^2))',0,n,n); % normalize columns
%   v = sqrt(0.1)*randn(m,1);
%   b = A*x0 + v;
%   fprintf('solving instance with %d examples, %d variables\n', m, n);
%   fprintf('nnz(x0) = %d; signal-to-noise ratio: %.2f\n', nnz(x0), norm(A*x0)^2/norm(v)^2);
%   gamma_max = norm(A'*b,'inf');
%   gamma = 0.01*gamma_max;
%   %Cached computations for all methods
%   AtA = A'*A;
%   Atb = A'*b;
L = max(eigs(AtA));
lambda = 1/L;
beta = 0.5;
%   %Global constants and defaults
%   MAX_ITER = 200;
%   ABSTOL  = eps;
%   
  f = @(u) 0.5*sum_square(double(A*u-b));

  % Run 
  [y0,x64] = entrypoint('double',f,x0,A,b,AtA,Atb,lambda,gamma,beta,MAX_ITER,ABSTOL);
  fprintf('Solved with double precision\n');
  [psnr_64,snr_64] =  psnr(x64,x0);
  ssimval_64 = ssim(x64,x0);
  
  [y01,x32] = entrypoint('single',f,x0,A,b,AtA,Atb,lambda,gamma,beta,MAX_ITER,ABSTOL);
  fprintf('Solved with single precision\n');
  [psnr_32,snr_32] =  psnr(x32,x0);
  ssimval_32 = ssim(x32,x0);
%   
  [y16,x16] = entrypoint('fixed16',f,x0,A,b,AtA,Atb,lambda,gamma,beta,MAX_ITER,ABSTOL); 
  fprintf('Solved with half precision\n');
  [psnr_16,snr_16] =  psnr(x16,x0);
  ssimval_16 = ssim(x16,x0);
% 
  [y12,x12] = entrypoint('fixed12',f,x0,A,b,AtA,Atb,lambda,gamma,beta,MAX_ITER,ABSTOL); 
  fprintf('Solved with 12 bits\n');
  [psnr_12,snr_12] =  psnr(x12,x0);
  ssimval_12 = ssim(x12,x0);  
%   
  [y8,x8] = entrypoint('fixed8',f,x0,A,b,AtA,Atb,lambda,gamma,beta,MAX_ITER,ABSTOL);
  fprintf('Solved with 8 bits\n');
  [psnr_8,snr_8] =  psnr(x8,x0);
  ssimval_8 = ssim(x8,x0);
% 
  [y2,x2] = entrypoint('fixed2',f,x0,A,b,AtA,Atb,lambda,gamma,beta,MAX_ITER,ABSTOL);
  fprintf('Solved with 2 bits\n');
  [psnr_2,snr_2] =  psnr(x2,x0);
  ssimval_2 = ssim(x2,x0);
% 
%   % Plot 
  figure;
  subplot(5,1,1);semilogy(y0,'k'); 
  legend('Ground truth objective') 
  title('Ground truth (Double precision)') 

  subplot(5,2,3);semilogy(y01,'k'); 
  title('Single precision floating-point output') 
  subplot(5,2,4);
  try
      semilogy(abs(y0(1:numel(y01))-double(y01)),'r');
      title('Residual error');
  catch
      semilogy(abs(y0-double(y01(1:numel(y0)))),'r');
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
%   
%   %pg_solv_mex
%   pg_solv_analyse(x0,A,b,AtA,Atb,lambda,gamma, beta, MAX_ITER, ABSTOL);
% 
  if (scenario == 'R')
      figure;
    %   imshow([reshape(x0,x0m,x0n),reshape(x64,x0m,x0n), ...
    %      reshape(x32,x0m,x0n),reshape(x16,x0m,x0n),reshape(x12,x0m,x0n),reshape(x8,x0m,x0n)],'InitialMagnification', 800);
      h64 = subplot(2,3,1), imshow(reshape(x64,x0m,x0n), 'Parent', h64), title(h64, ['BW = 64 bits (',num2str(round(100*64/64)),'%); SSIM = ', num2str(ssimval_64),'; PSNR = ', num2str(psnr_64)]); hold  on
      h32 = subplot(2,3,2), imshow(reshape(x32,x0m,x0n), 'Parent', h32), title(h32, ['BW = 32 bits (',num2str(round(100*32/64)),'%); SSIM = ', num2str(ssimval_32),'; PSNR = ', num2str(psnr_32)]); hold  on
      h16 = subplot(2,3,3), imshow(reshape(x16,x0m,x0n), 'Parent', h16), title(h16, ['BW = 16 bits (',num2str(round(100*16/64)),'%); SSIM = ', num2str(ssimval_16),'; PSNR = ', num2str(psnr_16)]); hold  on
      h12 = subplot(2,3,4), imshow(reshape(x12,x0m,x0n), 'Parent', h12), title(h12, ['BW = 12 bits (',num2str(round(100*12/64)),'%); SSIM = ', num2str(ssimval_12),'; PSNR = ', num2str(psnr_12)]); hold  on
      h8 = subplot(2,3,5), imshow(reshape(x8,x0m,x0n), 'Parent', h8), title(h8, ['BW = 8 bits (',num2str(round(100*8/64)),'%); SSIM = ', num2str(ssimval_8),'; PSNR = ', num2str(psnr_8)]); hold  on
      h2 = subplot(2,3,6), imshow(reshape(x2,x0m,x0n), 'Parent', h2), title(h2, ['BW = 2 bits (',num2str(round(100*2/64)),'%); SSIM = ', num2str(ssimval_2),'; PSNR = ', num2str(psnr_2)]); hold  off
  end
end 
