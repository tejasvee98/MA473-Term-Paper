function [U,x,t]=piecewise_spatial_mesh_bs_american(S_max,T,N,M,r,sigma,K,alpha,beta_star,epsilon,bc1,bc2,C,q,mu,method)
    x=zeros(1,N+1);
    h=zeros(1,N+1);
    for i=2:N+1
        x(i)=get_x(i-1, N, K, alpha, beta_star, S_max, epsilon);
        h(i)=get_h(i-1, N, K, alpha, beta_star, S_max, epsilon);
    end
    t=(0:T/M:T);
    tau=T/M;
    U=zeros(N+1,M+1);
    U(1,1:M)=bc1(t);
    U(N+1,1:M)=bc2(t);
    for i=1:N+1
        U(i,M+1)=pi_epsilon(K-x(i),epsilon);
    end
    for j=M:-1:1
        b1=U(1,j);
        b2=U(N+1,j);
        if method=="fsolve"
            F=@(v) BSP_Operator_with_Jacobian(v,U(:,j+1),N,b1,b2,tau,t(j),r,sigma,h,x,C,q,mu);
            options = optimoptions('fsolve');
            % doesn't reduce iterations Works without it too as well
            %options.SpecifyObjectiveGradient=true;
            %options.Display='iter';
            % 'trust region' works just as well
            options.Algorithm = 'levenberg-marquardt';
            U(:,j) = fsolve(F,U(:,j+1),options);
        elseif method == "newton-raphson"      
        % Newton raphson approach mentioned in paper.
            F=@(v) BSP_Operator(v,U(:,j+1),N,b1,b2,tau,t(j),r,sigma,h,x,C,q,mu);
            J=@(v) BSP_Operator_Jacobian(v,N,tau,t(j),r,sigma,h,x,C,q,mu);
            U(:,j) = my_newton_raphson(F,J,U(:,j+1),1e-5);
        end
        %U(1,j)=b1;
        %U(N+1,j)=b2;
    end
end