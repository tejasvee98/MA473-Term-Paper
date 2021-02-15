clc; clear all;
format long 
T=1; r=@(t) 0.06; sigma = @(x,t) 0.4*(2 + sin(x)); K=25; S_max=100;
M=64;N=64;
epsilon=1e-4;
% sigma^2 >= alpha > 0 and beta_star >= r >= beta > 0
alpha = 0.01;
beta_star = 100;
% add smooth tc
%tc=@(x) max(x-K,0);
tc=@(x) pi_epsilon(x-K,epsilon);
bc1=@(t) 0;
bc2=@(t) S_max - K .* exp(-r(0).*(T-t));
[U,x,t]=piecewise_spatial_mesh_bs(S_max,T,N,M,r,sigma,K,alpha,beta_star,epsilon,tc,bc1,bc2);

[X,T]=meshgrid(x,t);
surf(X,T,U');
xlabel('x');
ylabel('t');
zlabel('U(x,t)');
