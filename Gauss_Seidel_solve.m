function x=Gauss_Seidel_solve(A,b,x0,eps)
    D=diag(diag(A));
    L=tril(A)-D;
    U=triu(A)-D;
    x_new=x0;
    x_old=x0;
    D_inv=inv(D);
    while abs(max(A*x_new-b))>eps
        x_new=D_inv*(b-U*x_old-L*x_new);
        x_old=x_new;
    end
    x=x_new(:,1);
end