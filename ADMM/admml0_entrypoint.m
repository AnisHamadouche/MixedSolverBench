function [y,x,z,u,y_ref,x_ref,z_ref,u_ref] = admml0_entrypoint(A, b, Atb, lambda, gamma, rho, MAX_ITER, m, n, ABSTOL, RELTOL,dt)
    T = mytypes(dt);  
    %input
    A=cast(A, 'like', T.x);
    b=cast(b, 'like', T.x);
    AtA=cast(A'*A, 'like', T.x);
    Atb=cast(Atb, 'like', T.x);
    %parameters
    lambda=cast(lambda, 'like', T.param);
    rho = 1/lambda;
    gamma=cast(gamma, 'like', T.param);
    %beta=cast(beta, 'like', T.param);
    
    h = struct();
    h = admm_proxl0_solv(h,A, b, Atb, lambda, gamma, rho, MAX_ITER, m, n, ABSTOL, RELTOL,T);
    y = h.p_admm;
    x = double(h.x_admm);
    z = double(h.z_admm);
    u = double(h.u_admm);
    
    T = mytypes('double');
    h_ref = struct();
    h_ref = admm_proxsolv(h,A, b, Atb, lambda, gamma, rho, MAX_ITER, m, n, ABSTOL, RELTOL,T);
    y_ref = h_ref.p_admm;
    x_ref = double(h_ref.x_admm);
    z_ref = double(h_ref.z_admm);
    u_ref = double(h_ref.u_admm);
end