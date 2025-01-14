function [L, U] = factor(A, rho) % (H, rho, lambda_x, A, L_k, Mx)
    %A = double(A);
    [m, n] = size(A);
    if ( m >= n )    % if skinny
       L = chol( double(A'*A + rho*eye(n)), 'lower' ); % rho*eye(n) => WLM-ADMM => lambda_x*rho*A'*L_k*A+lambda_x*Mx
    else            % if fat
       L = chol( double(eye(m) + 1/rho*(A*A')), 'lower' );
    end

    % force matlab to recognize the upper / lower triangular structure
    L = L;
    U = L';
end
