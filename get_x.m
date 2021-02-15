function x_i = get_x(i, N, K, alpha, beta_star, S_max, epsilon)
    h = (K-epsilon)/(1 + alpha*(N/4 - 2)/beta_star);
    if i==0
        x_i = 0;
    elseif i==1
        x_i = h;
    elseif i < N/4
        x_i = h.*(1 + (i-1).*alpha/beta_star);
    elseif i==(N/4)
        x_i = K;
    elseif i == (N/4 + 1)
        x_i = K + epsilon;
    else
        x_i = K + epsilon + (S_max - K - epsilon).*(i - N/4 - 1)/(3*N/4 - 1);
    end   
end