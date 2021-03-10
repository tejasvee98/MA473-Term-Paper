function [U,x,t]=piecewise_spatial_mesh_bs(S_max,T,N,M,r,sigma,K,alpha,beta_star,epsilon,tc,bc1,bc2)
    x=zeros(1,N+1);
    h=zeros(1,N+1);
    for i=2:N+1
        x(i)=get_x(i-1, N, K, alpha, beta_star, S_max, epsilon);
        h(i)=get_h(i-1, N, K, alpha, beta_star, S_max, epsilon);
    end
    t=(0:T/M:T);
    tau=T/M;
    U=zeros(N+1,M+1);
    U(:,M+1)=tc(x);
    U(1,:)=bc1(t);
    U(N+1,:)=bc2(t);
    A=zeros(N+1,N+1);
    A(1,1)=1;
    A(1,2)=0;
    A(N+1,N)=0;
    A(N+1,N+1)=1;
    C=zeros(N+1,1);
    for j=M:-1:1
        for i=2:N
            x_i = x(i);
            h_i = h(i);
            h_i_plus_1 = h(i+1);
            sig_ij = sigma(x_i,t(j))^2;
            r_j = r(t(j));
            A(i,i-1) = (tau * x_i) * (r_j - (sig_ij * x_i / (h_i))) / (h_i + h_i_plus_1);
            A(i,i) = 1 + (sig_ij  * tau * (x_i^2))/(h_i * h_i_plus_1) + r_j * tau;
            A(i,i+1) = -((tau * x_i) * (r_j + (sig_ij * x_i / (h_i_plus_1))) / (h_i + h_i_plus_1));
        end
        % Can also use iterative methods to solve Ax=b
        C(1)=U(1,j);
        C(N+1)=U(N+1,j);
        C(2:N)=U(2:N,j+1);
        U(:,j)=A\C;
    end
end