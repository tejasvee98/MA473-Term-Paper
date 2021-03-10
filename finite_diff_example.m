clc; clear all;
format long 
% Use different values of Smax to get a better idea
T=1; rf=@(t) 0.06; sigma = @(x,t) 0.2*(1 + (t * exp(-x))); K=25; S_max=100;
% sigma = @(x,t) 0.4*(2+sin(x));
M=64;N=64;
h=S_max/M;
k=-T/N;

tc=@(x) max(x-K,0);
bc1=@(t) 0;
bc2=@(t) S_max - K .* exp(-rf(0).*(1-t));

p=@(x,t) 0.5*(sigma(x,t)^2);
q=@(x,t) -rf(t)*x;
r=@(x,t) -rf(t);
s=@(x,t) 0;
% use crank nicholson as well to check
[U,x,t,~,~]=btcs(h,k,0,S_max,T,0,p,q,r,s,tc,bc1,bc2);
[X,T]=meshgrid(x,t);
figure
surf(X,T,U');
xlabel('x');
ylabel('t');
zlabel('U(x,t)');
