function J=BSP_Operator_Jacobian(v,N,tau,t_j,r,sigma,h,x,C,q,mu)
    J=zeros(N+1,N+1);
    J(1,1)=1;
    J(N+1,N+1)=1;
    r_j = r(t_j);
    for i=2:N
        x_i = x(i);
        h_i = h(i);
        h_i_plus_1 = h(i+1);
        sig_ij = sigma(x_i,t_j)^2;
        c1=(tau * x_i) * (r_j - (sig_ij * x_i / (h_i))) / (h_i + h_i_plus_1);
        c2=1 + (sig_ij  * tau * (x_i^2))/(h_i * h_i_plus_1) + r_j * tau;
        c3=-((tau * x_i) * (r_j + (sig_ij * x_i / (h_i_plus_1))) / (h_i + h_i_plus_1));
        J(i,i-1)=c1;
        J(i,i)=c2 + (C*mu)/((v(i)+mu - q(x_i))^2);
        J(i,i+1)=c3;
    end
end