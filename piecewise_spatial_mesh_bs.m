function [U,s,t]=piecewise_spatial_mesh_bs(S_max,T,N,M,r,sigma,K,alpha,beta_star,epsilon,tc,bc1,bc2)
    s=(0:S_max/N:S_max);
    t=(0:T/M:T);
    U=zeros(N+1,M+1);
    U(:,M+1)=tc(s);
    U(1,:)=bc1(t);
    U(N+1,:)=bc2(t);
    A=zeros(N+1,N+1);
    A(1,1)=1;
    A(1,2)=0;
    A(N+1,N)=0;
    A(N+1,N+1)=1;
    for j=M:-1:1
        for i=2:N
            x_i = get_x(i, N, K, alpha, beta_star, S_max, epsilon);
            h_i = get_h(i, N, K, alpha, beta_star, S_max, epsilon);
            h_i_plus_1 = get_h(i+1, N, K, alpha, beta_star, S_max, epsilon);
            sig_ij = sigma(i,j)^2;
            r_j = r(j);
            A(i,i-1) = (t(j+1) * x_i) * (r_j - (sig_ij * x_i / (h_i))) / (h_i + h_i_plus_1);
            A(i,i) = 1 + (sig_ij  * t(j+1) * (x_i^2))/(h_i * h_i_plus_1) + r(j) * t(j+1);
            A(i,i+1) = -((t(j+1) * x_i) * (r_j + (sig_ij * x_i / (h_i_plus_1))) / (h_i + h_i_plus_1));
        end
        % Can also use iterative methods to solve Ax=b
        U(:,j)=A\U(:,j+1);
    end
end