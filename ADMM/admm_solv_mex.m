%Test script for instrumenting and executing code and showing proposed 
%fixed-point types and potential overflows for application variables

if exist('pg_solv_mex.mexw64', 'file')==2
  delete('pg_solv_mex.mexw64');
end

% TEST INPUT 
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

x0 = rand(n,1);
z0 = x0;
u0 = zeros(n,1);

gamma_max = norm(H'*b,'inf');
gamma = 0.1*gamma_max;

% cached computations for all methods
HtH = H'*H;
Htb = H'*b;

%% Global constants and defaults

MAX_ITER = 100;
ABSTOL   = 1e-4;
RELTOL   = 1e-2;

%% Extra parameters for ADMM
lambda = 1;
rho = 1/lambda;

%% Build 
buildInstrumentedMex admm_solv ... 
  -args {x0,z0,u0,H, b, Htb, lambda, gamma, rho, MAX_ITER, m, n, ABSTOL, RELTOL} -histogram 

%% Run 
y = admm_solv_mex(x0,z0,u0,H, b, Htb, lambda, gamma, rho, MAX_ITER, m, n, ABSTOL, RELTOL);

%% Show 
showInstrumentationResults admm_solv_mex ... 
  -defaultDT numerictype(1,16) -proposeFL 