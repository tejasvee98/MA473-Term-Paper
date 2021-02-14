clc; clear all;
T=1; r=@(t) 0.06; sigma = @(x,t) 0.2*(1 + (t * exp(-x))); K=25; S_max=100;
M=64;N=64;
epsilon=1e-4;
% sigma^2 >= alpha > 0 and beta_star >= r >= beta > 0
alpha = 0.01;
beta_star = 1;

tc=@(x) pi_epsilon(x-K,epsilon);
bc1=@(t) 0;
bc2=@(t) S_max - K * exp(-r(0));

[U,s,t]=piecewise_spatial_mesh_bs(S_max,T,N,M,r,sigma,K,alpha,beta_star,epsilon,tc,bc1,bc2);

[S,T]=meshgrid(s,t);
surf(S,T,U');
xlabel('x');
ylabel('t');
zlabel('U(x,t)');
