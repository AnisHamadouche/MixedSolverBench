function p = objective(A, b, gamma, x, z)
    p = 0.5*sum_square(A*x - b)+ gamma*l1_norm(z);
end
