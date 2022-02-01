%ADMM Solver
% lambda = 1;
% rho = 1/lambda;
% h = admm(h, A, b, Atb, lambda, gamma, rho, MAX_ITER, m, n, ABSTOL, RELTOL);
function y = admm_proxsolv(x0,z0,u0,A, b, Atb, lambda, gamma, rho, MAX_ITER, m, n, ABSTOL, RELTOL,dt)
    tic;
    %Data type casting
    T = mytypes(dt);
    x = cast(x0, 'like', T.x);%zeros(n,1, 'like', T.x);
    z = cast(z0, 'like', T.x);%zeros(n,1, 'like', T.x);
    zold = zeros(n,1, 'like', T.x);
    u = cast(u0, 'like', T.x);% zeros(n,1, 'like', T.x);
    q = zeros(n,n, 'like', T.x);
    [L, U] = factor(double(A),double(rho));
    L = cast(L,'like',T.x);
    U = cast(U,'like',T.x);
    A = cast(A,'like',T.x);
    for k = 1:MAX_ITER
        admm_p_i(:,k)   = objfunc(A, b, gamma, x, z);
        % x-update
        q = Atb + rho*(z - u);
        if m >= n
           x(:) = U \ (L \ q);
        else
           x(:) = lambda*(q - lambda*(double(A')*(double(U) \ ( double(L) \ (double(A*q)) ))));
        end
        % z-update
        zold(:) = z;
        z(:) = prox_l1(x + u, lambda*gamma);
        % u-update
        u(:) = u + x - z;
        % diagnostics, reporting, termination checks
        r_i(:,k)   = norm(double(x) - double(z));
        s_i(:,k)   = norm(-double(rho)*(double(z) - double(zold)));
        eps_pri(:,k)  = sqrt(n)*ABSTOL + RELTOL*max(norm(double(x)), norm(-double(z)));
        eps_dual(:,k) = sqrt(n)*ABSTOL + RELTOL*norm(double(rho)*double(u));
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