%Test script for instrumenting and executing code and showing proposed 
%fixed-point types and potential overflows for application variables

if exist('pg_solv_mex.mexw64', 'file')==2
  delete('pg_solv_mex.mexw64');
end

% TEST INPUT 
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

%% Build 
buildInstrumentedMex admm_solv ... 
  -args {x0,z0,u0,A, b, Atb, lambda, gamma, rho, MAX_ITER, m, n, ABSTOL, RELTOL} -histogram 

%% Run 
y = admm_solv_mex(x0,z0,u0,A, b, Atb, lambda, gamma, rho, MAX_ITER, m, n, ABSTOL, RELTOL);

%% Show 
showInstrumentationResults admm_solv_mex ... 
  -defaultDT numerictype(1,16) -proposeFL 