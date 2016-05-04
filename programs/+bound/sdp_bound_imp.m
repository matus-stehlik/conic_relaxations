function [lb, ub, xf] = sdp_bound_imp(W,P)
% returns lower bound lb, upper bound ub and the feasible solution xu
% which generates ub for the problem
% min x^TWx, s.t. x is in {-1,1}^n, 
% and some coordinates of x, x(P.neg) = -1; x(P.pos) = 1 are given


% instead of solving min xWx with constraints on coeffs of x,
% we will solve problem with these values already substituted
% reducing the dimension 

N = size(W,1);
xk = [ones(length(P.pos),1);-ones(length(P.neg),1)];
K = [P.pos,P.neg];  % known coordinates of x
U = setdiff(1:N,K); % unknown coords of x
n = length(U)+1;    % dimension of variable Y in sdp relax
w0 = xk'*W(K,K)*xk; % M(1,1)
wk = W(U,K)*xk;     % M(2:end,1)

Q = [w0,wk';wk,W(U,U)];
lambda = min(min(abs(Q(Q>0))))/(n-1);

% if all the variables are set, we dont need any optimization
if numel(K) == N,
    lb = w0;
    ub = w0;
    xf = xk;
    return;
end


% min x^TWx with constraints on x given in P is equivalent to 
% min xu^TW(U,U)xu + 2wk^Tx + w0, s.t. xu is from {-1,1}^(n-1)
% which can be simplified to min trace(MY), s.t. diag(Y)=1 and Y is PSD
% where M = [w0,wk';wk,W(U,U)]

cvx_begin sdp
    cvx_quiet(true);
    variable Y(n,n) symmetric
    variable z(1,1)
    minimize ( trace( [w0,wk';wk,W(U,U)]*Y) + lambda*z) 
    diag(Y) == 1
    Y*ones(n,1) >= -z
    Y*ones(n,1) <= z
    Y == semidefinite(n);
cvx_end

z

lb = ceil(cvx_optval)  ;  % from sdp we have obtained lower bound
lb = cvx_optval  ; 
                            % !!! this should be done in more clever way
% xu = bound.triv_ubound(Y(2:end,2:end));
xu = bound.GW_ubound(Y,W(U,U),100);

xf = zeros(N,1);
xf(U) = xu;
xf(K) = xk;
ub1 = xu'*W(U,U)*xu + 2*wk'*xu+ w0 ;
ub2 = ub1 - 4*wk'*xu;
[ub,ind] = min([ub1,ub2]);
if ind==2, xf(U) = -xf(U); end

end