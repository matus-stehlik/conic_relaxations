function [optVal] = sdp(A,n)

cvx_begin sdp
    cvx_quiet(true);
    variable X(n,n) symmetric
    minimize (trace(A'*X))
    diag(X) == 1
    X == semidefinite(n);
cvx_end

optVal = cvx_optval;

end