function admm_benchmarking
clear;
clc;
close all;
addpath 'C:\Users\ah225\Documents\MATLAB\RPFOSolver Benchmark\ADMM'

%% Problem data

% s = RandStream.create('mt19937ar','seed',0);
% RandStream.setDefaultStream(s);

m = 100;       % number of examples
n = 500;      % number of features

x0 = randn(n,1);
z0 = randn(n,1);
u0 = randn(n,1);
A = randn(m,n);
A = A*spdiags(1./sqrt(sum(A.^2))',0,n,n); % normalize columns
v = sqrt(0.001)*randn(m,1);
b = A*x0 + v;

fprintf('solving instance with %d examples, %d variables\n', m, n);
fprintf('nnz(x0) = %d; signal-to-noise ratio: %.2f\n', nnz(x0), norm(A*x0)^2/norm(v)^2);

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
  y0 = admm_entrypoint(x0,z0,u0,A, b, Atb, lambda, gamma, MAX_ITER, m, n, ABSTOL, RELTOL,'double');
  y01 = admm_entrypoint(x0,z0,u0,A, b, Atb, lambda, gamma, MAX_ITER, m, n, ABSTOL, RELTOL,'single');
  y8 = admm_entrypoint(x0,z0,u0,A, b, Atb, lambda, gamma, MAX_ITER, m, n, ABSTOL, RELTOL,'fixed8');
  y16 = admm_entrypoint(x0,z0,u0,A, b, Atb, lambda, gamma, MAX_ITER, m, n, ABSTOL, RELTOL,'fixed16');
  y12 = admm_entrypoint(x0,z0,u0,A, b, Atb, lambda, gamma, MAX_ITER, m, n, ABSTOL, RELTOL,'fixed12');

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

  admm_solv_mex
end 