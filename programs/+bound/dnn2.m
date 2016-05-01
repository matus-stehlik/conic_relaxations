function [lb, ub, xf, utime, ltime] = dnn2(W,P)
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
n = length(U);    % dimension of variable Y in sdp relax
w0 = xk'*W(K,K)*xk; % M(1,1)
wk = W(U,K)*xk;     % M(2:end,1)

% if all the variables are set, we dont need any optimization
if numel(K) == N,
    utime = -1;
    ltime = -1;
    ub = w0;
    lb = w0;
    xf = xk;
    return;
end


% min x^TWx with constraints on x given in P is equivalent to 
% min xu^TW(U,U)xu + 2wk^Tx + w0, s.t. xu is from {-1,1}^(n-1)
% which can be simplified to min trace(MY), s.t. diag(Y)=1 and Y is PSD
% where M = [w0,wk';wk,W(U,U)]
cvx_solver sedumi;
tic;
cvx_begin sdp
    cvx_quiet(true);
    variable Y(n,n) symmetric
    variable UU(n,n) symmetric
    variable VV(n,n)
    variable y(n,1) 
    variable u(n,1)
    
    maximize ( 4*w0 + 8*wk'*y + 4*trace(W(U,U)*Y) ) 
    diag(Y) == y
    diag(UU) == u
    diag(VV) == 0
    y + u == 1
    y >= 0
    u >= 0
    Y(:) >= 0
    UU(:) >= 0
    VV(:) >= 0
    [1 y' u'; y Y VV';u VV UU]  == semidefinite(2*n+1);
cvx_end

ub = cvx_optval  ;  % from sdp we have obtained upper bound
utime = toc;
tic;
xl = bound.GW_lbound(2*[1 y'; y Y]+1,W(U,U),100);
% xl = bound.GW_lbound(4*[1 y'; y Y] - 2*(ones(n+1,1)*[1 Y(1,:)] + ...
%                     [1;Y(:,1)]*ones(1,n+1)) + 1,W(U,U),100);

xf = zeros(N,1);
xf(U) = xl;
xf(K) = xk;
lb1 = xl'*W(U,U)*xl + 2*wk'*xl+ w0 ;
lb2 = lb1 - 4*wk'*xl;
[lb,ind] = max([lb1,lb2]);
if ind==2, xf(U) = -xf(U); end
ltime = toc;
end