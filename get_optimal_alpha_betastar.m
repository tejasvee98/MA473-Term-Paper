% sigma^2 >= alpha > 0 and beta_star >= r >= beta > 0
function [alpha, beta_star]=get_optimal_alpha_betastar(S_max,T,N,M,r,sigma,K,epsilon,tc,bc1,bc2)
    alpha = 0.01;
    beta_star = 1;
    min_error = 1e5;
    for a=0.01:0.001:0.04
        [U_exact,x_exact,t_exact]=piecewise_spatial_mesh_bs(S_max,T,2048,M,r,sigma,K,a,beta_star,epsilon,tc,bc1,bc2);
        [X,Y]=meshgrid(x_exact,t_exact);
        V=@(x,y) interp2(X,Y,U_exact',x,y);
        [U,x,t]=piecewise_spatial_mesh_bs(S_max,T,N,M,r,sigma,K,a,beta_star,epsilon,tc,bc1,bc2);
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
        [U_exact,x_exact,t_exact]=piecewise_spatial_mesh_bs(S_max,T,2048,M,r,sigma,K,alpha,b,epsilon,tc,bc1,bc2);
        [X,Y]=meshgrid(x_exact,t_exact);
        V=@(x,y) interp2(X,Y,U_exact',x,y);
        [U,x,t]=piecewise_spatial_mesh_bs(S_max,T,N,M,r,sigma,K,alpha,b,epsilon,tc,bc1,bc2);
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