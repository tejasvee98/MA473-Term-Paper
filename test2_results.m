clc;
format long e
N=[128,256,512,1024]; M=1024;

T=1; r=@(t) 0.06; sigma = @(x,t) 0.4*(2+sin(x)); K=25; S_max=100;
epsilon=1e-4;

alpha = 0.03;
beta_star = 0.1;

bc1=@(t) 0;
bc2=@(t) S_max - K .* exp(-r(0).*(1-t));

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