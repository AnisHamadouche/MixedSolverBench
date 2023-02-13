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
    m = 100;
    noise_mean = 0;
    noise_var = 0.0001; % 10^{-5}
    inputImage = im2double(imread('./data/einstein.jpg'));
    % reduce resolution to fit in memory
    N = 10;
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
    %imshow([reshape(x0,x0m,x0n),reshape(b,x0m,x0n)]);
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
  [y01,x32] = entrypoint('single',f,x0,A,b,AtA,Atb,lambda,gamma,beta,MAX_ITER,ABSTOL);
  fprintf('Solved with single precision\n');
  [y16,x16] = entrypoint('fixed16',f,x0,A,b,AtA,Atb,lambda,gamma,beta,MAX_ITER,ABSTOL); 
  fprintf('Solved with half precision\n');
  [y12,x12] = entrypoint('fixed12',f,x0,A,b,AtA,Atb,lambda,gamma,beta,MAX_ITER,ABSTOL); 
  fprintf('Solved with 12 bits\n');
  [y8,x8] = entrypoint('fixed8',f,x0,A,b,AtA,Atb,lambda,gamma,beta,MAX_ITER,ABSTOL);
  fprintf('Solved with 8 bits\n');
  % Plot 
  subplot(4,1,1);plot(y0,'k'); 
  legend('Ground truth objective') 
  title('Ground truth') 

  subplot(4,2,3);plot(y01,'k'); 
  title('Single precision floating-point output') 
  subplot(4,2,4);
  try
      semilogy(abs(y0(1:numel(y01))-double(y01)),'r');
  catch
      semilogy(abs(y0-double(y01(1:numel(y0)))),'r');
  end
  title('Single precision floating-point error') 

  subplot(4,2,5);plot(y12,'k'); 
  title('12-bit fixed-point output') 
  subplot(4,2,6);
  try
      semilogy(abs(y0(1:numel(y12))-double(y12)),'r');
  catch
      semilogy(abs(y0-double(y12(1:numel(y0)))),'r');
  end
  title('12-bit fixed-point error') 

  subplot(4,2,7);plot(y16,'k'); 
  title('16-bit fixed-point output') 
  xlabel('Iterations, k') 
  subplot(4,2,8);
  try
      semilogy(abs(y0(1:numel(y16))-double(y16)),'r');
  catch
      semilogy(abs(y0-double(y16(1:numel(y0)))),'r');
  end
  title('16-bit fixed-point error') 
  xlabel('Iterations, k') 
  
  %pg_solv_mex
  pg_solv_analyse(x0,A,b,AtA,Atb,lambda,gamma, beta, MAX_ITER, ABSTOL);

  figure;
  imshow([reshape(x0,x0m,x0n),reshape(x64,x0m,x0n), ...
      reshape(x32,x0m,x0n),reshape(x16,x0m,x0n),reshape(x12,x0m,x0n),reshape(x8,x0m,x0n)],'InitialMagnification', 800);
end 