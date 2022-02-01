%% TEST INPUT  
f = @(u) 0.5*sum_square(double(A*u-b));
m = 50;       % number of examples
n = 10;      % number of features
%x0 = sprandn(n,1,0.05);
x0 = randn(n,1);
A = randn(m,n);
%A = A*spdiags(1./sqrt(sum(A.^2))',0,n,n); % normalize columns
v = sqrt(0.01)*randn(m,1);
b = A*x0 + v;
fprintf('solving instance with %d examples, %d variables\n', m, n);
fprintf('nnz(x0) = %d; signal-to-noise ratio: %.2f\n', nnz(x0), norm(A*x0)^2/norm(v)^2);
gamma_max = norm(A'*b,'inf');
gamma = 0.01*gamma_max;
%% Cached computations for all methods
AtA = A'*A;
Atb = A'*b;
L = max(eigs(AtA));
lambda = 1/L;
beta = 0.5;

%% Global constants and defaults
MAX_ITER = 200;
ABSTOL   = eps;

% ALGORITHM 
y = pg_solv(f,x0,A,b,AtA,Atb,lambda,gamma,beta,MAX_ITER,ABSTOL); 

% VERIFY RESULTS 
yExpected = pg_proxsolv(f,x0,A,b,AtA,Atb,lambda,gamma,beta,MAX_ITER,ABSTOL);
try
    plot(y.p_i-yExpected.p_i(1:numel(y.p_i))) 
catch
    plot(y.p_i(1:numel(yExpected.p_i))-yExpected.p_i)
end
title('Error')