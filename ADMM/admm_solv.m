%ADMM Solver
% lambda = 1;
% rho = 1/lambda;
% h = admm(h, A, b, Atb, lambda, gamma, rho, MAX_ITER, m, n, ABSTOL, RELTOL);
function y = admm_solv(x0,z0,u0,A, b, Atb, lambda, gamma, rho, MAX_ITER, m, n, ABSTOL, RELTOL)
    tic;
    x = x0;
    z = z0;
    u = u0;
    [L, U] = factor(A, rho);

    admm_p_i   = zeros(1,MAX_ITER);
    r_i   = zeros(1,MAX_ITER);
    s_i   = zeros(1,MAX_ITER);
    eps_pri  = zeros(1,MAX_ITER);
    eps_dual = zeros(1,MAX_ITER);
    for k = 1:MAX_ITER
        % x-update
        q = Atb + rho*(z - u);
%         if m >= n
%            x = U \ (L \ q);
%         else
        x = lambda*(q - lambda*(A'*(U \ ( L \ (A*q) ))));
%         end
        % z-update
        zold = z;
        z = prox_l1(x + u, lambda*gamma);
        % u-update
        u = u + x - z;
        % diagnostics, reporting, termination checks
        admm_p_i(:,k)   = objfunc(A, b, gamma, x, z);
        r_i(:,k)   = norm(x - z);
        s_i(:,k)   = norm(-rho*(z - zold));
        eps_pri(:,k)  = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
        eps_dual(:,k) = sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
        if r_i(:,k) < eps_pri(:,k) && s_i(:,k) < eps_dual(:,k)
            admm_k = k; 
            break;
        end
    end
    x_admm = z;
    p_admm = admm_p_i(end);
    admm_time = toc;
    y = admm_p_i;
end