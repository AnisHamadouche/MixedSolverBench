function y = admm_entrypoint(x0,z0,u0,A, b, Atb, lambda, gamma, MAX_ITER, m, n, ABSTOL, RELTOL,dt)
    T = mytypes(dt);  
    %input
    A=cast(A, 'like', T.x);
    b=cast(b, 'like', T.x);
    %AtA=cast(AtA, 'like', T.x);
    Atb=cast(Atb, 'like', T.x);
    %parameters
    lambda=cast(lambda, 'like', T.param);
    rho = 1/lambda;
    gamma=cast(gamma, 'like', T.param);
    %beta=cast(beta, 'like', T.param);
    h = admm_proxsolv(x0,z0,u0,A, b, Atb, lambda, gamma, rho, MAX_ITER, m, n, ABSTOL, RELTOL,dt);
    y = h;
end