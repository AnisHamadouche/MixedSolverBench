function h = pg_proxsolv(f,x0,A,b,AtA,Atb,lambda, ...
                          gamma, beta, MAX_ITER,ABSTOL,dt)
        tic;
%       %Data type casting
        T = mytypes(dt); 
%         x = cast(x0, 'like', T.x);
%         x0 = cast(x0, 'like', T.x);
%         %output
%         h.x_i(:,1) = zeros(size(x0), 'like', T.x);
%         h.p_i(1) = zeros(1, 'like', T.y);
%         %input
%         A=cast(A, 'like', T.x);
%         b=cast(b, 'like', T.x);
%         AtA=cast(AtA, 'like', T.x);
%         Atb=cast(Atb, 'like', T.x);
%         %parameters
%         lambda=cast(lambda, 'like', T.param);
%         gamma=cast(gamma, 'like', T.param);
%         beta=cast(beta, 'like', T.param);
%         %internal
        z=zeros(size(x0), 'like', T.x);
        grad_x=zeros(size(x0), 'like', T.x);
        x = cast(x0, 'like', T.x);
        for k = 1:MAX_ITER
            grad_x(:) = AtA*x - Atb;
            z(:) = prox_l1(x - lambda*grad_x, lambda*gamma);
%         while 1
%             grad_x(:) = AtA*x - Atb;
%             z(:) = prox_l1(x - lambda*grad_x, lambda*gamma);
%             if f(z) <= f(x) + dot(grad_x,(z - x)) + (1/(2*lambda))*sum_square(z - x)
%                 break;
%             end
%             lambda = beta*lambda;
%         end
            x = z;
            h.x_i(:,k) = x;

            h.p_i(:,k) = objfunc(A, b, gamma, x, x);
    
            if k > 1 && abs(h.p_i(k) ...
                            - h.p_i(k-1)) < ABSTOL
                break;
            end
        end
        h.k = k;
        h.time = toc;
end