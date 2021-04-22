clc; clear all;
format long 
% Use different values of Smax to get a better idea
T=1; r=@(t) 0.06; sigma = @(x,t) 0.2; K=25; S_max=100;
%sigma = @(x,t) 0.4*(2+sin(x));
M=64;N=64;
epsilon=1e-4;
% sigma^2 >= alpha > 0 and beta_star >= r >= beta > 0
alpha = 0.03;
beta_star = 0.1;
% added smooth tc
%tc=@(x) max(x-K,0);
bc1=@(t) K;
bc2=@(t) 0;
mu=1e-4;
C=r(1)*K;
q=@(x) K-x;
[U,x,t]=piecewise_spatial_mesh_bs_american(S_max,T,N,M,r,sigma,K,alpha,beta_star,epsilon,bc1,bc2,C,q,mu,"newton-raphson");
v_star= @(x) max(K-x,0);
figure;
U_max=v_star(x);
[X,T]=meshgrid(x,t);
surf(X,T,U');
xlabel('x');
ylabel('t');
zlabel('U(x,t)');
figure
surf(X,T,U'-U_max);
xlabel('x');
ylabel('t');
zlabel('V-V*');
figure
s_f=zeros(M+1,1);
for j=1:M+1
    s_f(j) = x(find(abs(U(:,j)-K*ones(N+1,1)+x')< 1e-2, 1, 'last'));
end
plot(t,s_f);