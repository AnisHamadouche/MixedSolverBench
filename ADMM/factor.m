function [L, U] = factor(A, rho)
    %A = double(A);
    [m, n] = size(A);
    if ( m >= n )    % if skinny
       L = chol( A'*A + rho*eye(n), 'lower' );
    else            % if fat
       L = chol( eye(m) + 1/rho*(A*A'), 'lower' );
    end

    % force matlab to recognize the upper / lower triangular structure
    L = L;
    U = L';
end