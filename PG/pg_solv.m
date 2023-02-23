function y = pg_solv(x0,A,b,AtA,Atb,lambda, ...
                                  gamma, beta, MAX_ITER, ABSTOL)
        x = zeros(size(x0));
        z = zeros(size(x0));
        grad_x = zeros(size(x0));
        x_i=zeros(size(x0,1),MAX_ITER);
        p_i=zeros(1,MAX_ITER);
        k=1;
        time=0;
        tic;
        x = x0;
        f = @(u) 0.5*sum_square(double(A*u-b));
        y = objfunc(A, b, gamma, x0, x0);
        for k = 1:MAX_ITER
            grad_x = AtA*x - Atb;
            while 1
                z = prox_l1(x - lambda*grad_x, lambda*gamma);
                
                % back-tracking
                if f(z) <= f(x) + dot(grad_x',(z - x)) ...
                           + (1/(2*lambda))*sum_square(z - x)
                    break;
                end
                lambda = beta*lambda;
            end
            x = z;
            x_i(:,k) = x;
            p_i(:,k) = objfunc(A, b, gamma, x, x);
            time = toc;
            if k > 1 && abs(p_i(k) ...
                            - p_i(k-1)) < ABSTOL
                %h.k = k
                y = p_i(k);
                break;
            end
        end
end

function p = objfunc(A, b, gamma, x, z)
    p = 0.5*sum_square(A*x - b) + gamma*l1_norm(z);
end