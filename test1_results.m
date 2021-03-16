clc;
format long e

fprintf('Running test 1\n');

N=[128,256,512,1024]; M=1024;

T=1; r=@(t) 0.06; sigma = @(x,t) 0.2*(1 + (t * exp(-x))); K=25; S_max=100;
epsilon=1e-4;

bc1=@(t) 0;
bc2=@(t) S_max - K .* exp(-r(0).*(1-t));
 
%[alpha, beta_star]=get_optimal_alpha_betastar(S_max,T,64,M,r,sigma,K,epsilon,bc1,bc2);
%the optimizer takes a long time to run so we ran it once and used the
%values obtained below
alpha = 0.031;
beta_star = 1;

error=zeros(length(N),1);
[U_exact,x_exact,t_exact]=piecewise_spatial_mesh_bs(S_max,T,2048,M,r,sigma,K,alpha,beta_star,epsilon,bc1,bc2);

[X,Y]=meshgrid(x_exact,t_exact);
V=@(x,y) interp2(X,Y,U_exact',x,y);

for i=1:length(N)
    [U,x,t]=piecewise_spatial_mesh_bs(S_max,T,N(i),M,r,sigma,K,alpha,beta_star,epsilon,bc1,bc2);
    U=U';
    [x,t]=meshgrid(x,t);
    error(i)=max(abs(V(x,t)-U),[],'all');
end
R=V(x,t);
error
rate=zeros(length(N)-1,1);
for i=1:length(rate)
    rate(i)=log2(error(i)/error(i+1));
end
rate

% sigma^2 >= alpha > 0 and beta_star >= r >= beta > 0
function [alpha, beta_star]=get_optimal_alpha_betastar(S_max,T,N,M,r,sigma,K,epsilon,bc1,bc2)
    alpha = 0.01;
    beta_star = 1;
    min_error = 1e5;
    for a=0.01:0.001:0.04
        [U_exact,x_exact,t_exact]=piecewise_spatial_mesh_bs(S_max,T,2048,M,r,sigma,K,a,beta_star,epsilon,bc1,bc2);
        [X,Y]=meshgrid(x_exact,t_exact);
        V=@(x,y) interp2(X,Y,U_exact',x,y);
        [U,x,t]=piecewise_spatial_mesh_bs(S_max,T,N,M,r,sigma,K,a,beta_star,epsilon,bc1,bc2);
        U=U';
        [x,t]=meshgrid(x,t);
        error=max(abs(V(x,t)-U),[],'all');
        if error < min_error
            alpha = a;
            min_error = error;
        end
    end
    min_error = 1e5;
    for b=0.1:0.1:10
        [U_exact,x_exact,t_exact]=piecewise_spatial_mesh_bs(S_max,T,2048,M,r,sigma,K,alpha,b,epsilon,bc1,bc2);
        [X,Y]=meshgrid(x_exact,t_exact);
        V=@(x,y) interp2(X,Y,U_exact',x,y);
        [U,x,t]=piecewise_spatial_mesh_bs(S_max,T,N,M,r,sigma,K,alpha,b,epsilon,bc1,bc2);
        U=U';
        [x,t]=meshgrid(x,t);
        error=max(abs(V(x,t)-U),[],'all');
        if error < min_error
            beta_star = b;
            min_error = error;
        end
    end
    alpha
    beta_star
end