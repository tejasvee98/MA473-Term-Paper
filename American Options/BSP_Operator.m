function F=BSP_Operator(v,v_up,N,b1,b2,tau,t_j,r,sigma,h,x,C,q,mu)
    F=zeros(N+1,1);
    F(1)=v(1)-b1;
    F(N+1)=v(N+1)-b2;
    r_j = r(t_j);
    for i=2:N
        x_i = x(i);
        h_i = h(i);
        h_i_plus_1 = h(i+1);
        sig_ij = sigma(x_i,t_j)^2;
        c1=(tau * x_i) * (r_j - (sig_ij * x_i / (h_i))) / (h_i + h_i_plus_1);
        c2=1 + (sig_ij  * tau * (x_i^2))/(h_i * h_i_plus_1) + r_j * tau;
        c3=-((tau * x_i) * (r_j + (sig_ij * x_i / (h_i_plus_1))) / (h_i + h_i_plus_1));
        F(i)=c1*v(i-1) + c2*v(i) + c3*v(i+1) - (C*mu)/(v(i)+mu - q(x_i)) - v_up(i);
    end
end