%The proximal operator of the l1 norm with parameter lambda.
function x = prox_l1(v, lambda)
    x = max(0, v - lambda) - max(0, -v - lambda);
end 