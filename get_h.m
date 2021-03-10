function h_i = get_h(i, N, K, alpha, beta_star, S_max, epsilon)
    h = (K-epsilon)/(1 + alpha*(N/4 - 2)/beta_star);
    if i==1
        h_i = h;
    elseif i < N/4
        h_i = h*alpha/beta_star;
    elseif i < (N/4 + 2)
        h_i = epsilon;
    else
        h_i = (S_max - K - epsilon)/(3*N/4 - 1);
    end   
end