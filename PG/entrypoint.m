function y = entrypoint(dt,f,x0,A,b,AtA,Atb,lambda,gamma,...
                            beta,MAX_ITER,ABSTOL)
    T = mytypes(dt);  
    x0 = cast(x0, 'like', T.x);
    %input
    A=cast(A, 'like', T.x);
    b=cast(b, 'like', T.x);
    AtA=cast(AtA, 'like', T.x);
    Atb=cast(Atb, 'like', T.x);
    %parameters
    lambda=cast(lambda, 'like', T.param);
    gamma=cast(gamma, 'like', T.param);
    beta=cast(beta, 'like', T.param);
    h = pg_proxsolv(f,x0,A,b,AtA,Atb,lambda,gamma,beta,MAX_ITER,ABSTOL,dt);
    y=h.p_i;
end