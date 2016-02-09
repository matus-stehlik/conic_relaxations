function optVal = socp(A,n)

cvx_begin
    %cvx_quiet(true);

cvx_end

optVal = cvx_optval;

end