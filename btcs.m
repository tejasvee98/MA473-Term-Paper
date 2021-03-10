function [u,x,t,m,n]=btcs(h,k,x_a,x_b,t_a,t_b,p,q,r,s,ic,bc1,bc2)
    x=(x_a:h:x_b);
    t=(t_a:k:t_b);
    lam1=k/(h^2);
    lam2=k/(2*h);
    m=length(x);
    n=length(t);
    u=zeros(m,n);
    u(:,1)=ic(x);
    u(1,:)=bc1(t);
    u(m,:)=bc2(t);
    A=zeros(m,m);
    A(1,1)=1;
    A(1,2)=0;
    A(m,m-1)=0;
    A(m,m)=1;
    C=zeros(m,1);
    for j=2:n
        for i=2:m-1
            A(i,i-1)=p(x(i),t(j))*lam1 - q(x(i),t(j))*lam2;
            A(i,i)=1 - 2*p(x(i),t(j))*lam1 + k*r(x(i),t(j));
            A(i,i+1)=p(x(i),t(j))*lam1 + q(x(i),t(j))*lam2;
        end
       C(1) = u(1,j);
       C(m) = u(m,j);
       C(2:m-1)=u(2:m-1,j-1) + k*s(x(2:m-1),t(j));
       u(:,j)=A\C;
    end
end