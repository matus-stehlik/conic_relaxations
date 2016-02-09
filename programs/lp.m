function optVal = lp(A,n),
cvx_begin 
    cvx_quiet(true);
    variable X(n,n) symmetric
    variable t
    minimize (t)
    trace(A'*X)<= t
    diag(X) == 1
    vec(X) >= -1
    vec(X) <= 1;
    %X == semidefinite(n);

cvx_end

optVal = cvx_optval;
end