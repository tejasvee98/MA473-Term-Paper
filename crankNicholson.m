function [u,x,t,m,n]=crankNicholson(h,k,x_a,x_b,t_a,t_b,p,q,r,s,ic,bc1,bc2)
    % u_t + p*u_xx + q*u_x+r*u=s
    x=(x_a:h:x_b); t=(t_a:k:t_b);
    lam1=k/(h^2); lam2=(k/2*h);
    m=length(x); n=length(t);
    u=zeros(m,n);
    u(:,1)=ic(x);
    u(1,:)=bc1(t); u(m,:)=bc2(t);
    A=zeros(m,m); B=zeros(m,m); C=zeros(m,1);
    A(1,1)=1; A(1,2)=0;
    A(m,m-1)=0; A(m,m)=1;
    B(1,1)=1; B(1,2)=0;
    B(m,m-1)=0; B(m,m)=1;
    for j=2:n
        for i=2:m-1
            A(i,i-1)=0.5*(p(x(i),t(j))*lam1 - q(x(i),t(j))*lam2);
            A(i,i)=1 + 0.5*(-2*p(x(i),t(j))*lam1 + k*r(x(i),t(j)));
            A(i,i+1)=0.5*(p(x(i),t(j))*lam1 + q(x(i),t(j))*lam2);
        end
        for i=2:m-1
            B(i,i-1)=0.5*(-p(x(i),t(j-1))*lam1-q(x(i),t(j-1))*lam2);
            B(i,i)=1+0.5*(2*p(x(i),t(j-1))*lam1-k*r(x(i),t(j-1)));
            B(i,i+1)=0.5*(-p(x(i),t(j-1))*lam1+q(x(i),t(j-1))*lam2);
        end
        C=B*u(:,j-1) + 0.5*(k*s(x(i),t(j))+k*s(x(i),t(j-1)));
        C(1) = u(1,j); C(m) = u(m,j);
        u(:,j)=A\C;
    end
end