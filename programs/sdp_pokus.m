cvx_begin 
    variable X(n,n) symmetric
    variable t
    minimize (t)
    trace(A'*X)<= t
    diag(X) == 1
    vec(X) >= -1
    vec(X) <= 1;
    %X == semidefinite(n);

cvx_end