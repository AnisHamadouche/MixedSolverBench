function h = admm_proxsolv(h, A, b, Atb, lambda, gamma, rho, MAX_ITER, m, n, ABSTOL, RELTOL,T)

    tic;
    x = zeros(n,1, 'like', T.x);
    x_iter(:,1)=x;

    z = zeros(n,1, 'like', T.x);
    zold = zeros(n,1, 'like', T.x);
    z_iter(:,1)=z;

    u = zeros(n,1, 'like', T.x);
    u_iter(:,1)=u;

    q = zeros(n,1, 'like', T.x);

    [L, U] = factor(A, rho); % factor (A'A,rho)
    U_inv = pinv(U);
    L_inv = pinv(L);
    for k = 1:MAX_ITER
        % x-update
        q(:) = Atb + rho*(z - u);
        if m >= n
           %x = U \ (L \ q);
           x(:) = U_inv * (L_inv * q);
           x_iter(:,k)=x;
        else
           %x = lambda*(q - lambda*(A'*(U \ ( L \ (A*q) ))));
           x(:) = q/rho - (A'*(U_inv*(L_inv*(A*q))))/rho^2;
           x_iter(:,k)=x;
        end
        % z-update
        zold(:) = z;
        z(:) = prox_l1(x + u, lambda*gamma);
        z_iter(:,k)=z;
        
        % u-update
        u(:) = u + x - z;
        u_iter(:,k)=u;

        % diagnostics, reporting, termination checks
        h.admm_optval(k)   = objective(A, b, gamma, x, z);
        h.r_norm(k)   = norm(double(x - z));
        h.s_norm(k)   = norm(double(-rho*(z - zold)));
        h.eps_pri(k)  = sqrt(n)*ABSTOL + RELTOL*max(norm(double(x)), norm(double(-z)));
        h.eps_dual(k) = sqrt(n)*ABSTOL + RELTOL*norm(double(rho*u));
%         if h.r_norm(k) < h.eps_pri(k) && h.s_norm(k) < h.eps_dual(k)
%              break;
%         end
    end
    h.x_admm = x_iter;
    h.z_admm = z_iter;
    h.u_admm = u_iter;
    h.p_admm = h.admm_optval;
    h.admm_toc = toc;    
end
