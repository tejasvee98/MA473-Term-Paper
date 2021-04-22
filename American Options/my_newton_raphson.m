function x=my_newton_raphson(F,J,x0,tol)
    x=x0;
    n=1;
    x_new = x - J(x)\F(x);
    while(max(abs(x_new-x),[],'all')>tol)
        n=n+1;
        x=x_new;
        x_new = x - J(x)\F(x);
    end
end