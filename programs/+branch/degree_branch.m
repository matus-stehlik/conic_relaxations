function i = degree_branch(P,W),

N = size(W,1);
xk = [ones(length(P.pos),1);-ones(length(P.neg),1)];
K = [P.pos,P.neg];  % known coordinates of x
if numel(K) == N, 
    i = -1;
    return;
end

U = setdiff(1:N,K); % unknown coords of x
[~,ind] = max(diag(W(U,U)));
i = U(ind);


end