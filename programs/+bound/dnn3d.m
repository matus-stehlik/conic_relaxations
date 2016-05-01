function [lb, ub, xf, utime, ltime] = dnn3d(W,P,lam)
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



H1 = [-ones(n,1) eye(n) eye(n)]'*[-ones(n,1) eye(n) eye(n)] ...
    + [zeros(1,2*n+1); zeros(n,n+1) eye(n); zeros(n,1) eye(n) zeros(n)];

Q0  = 4*[w0 wk' zeros(1,n); wk W(U,U) zeros(n); zeros(n,2*n+1)];
H0 = sparse(1,1,1,2*n+1,2*n+1);
cvx_solver sedumi;

tic;
cvx_begin sdp
    cvx_quiet(true);
    variable y(1,1)
    variable NN(2*n+1,2*n+1) symmetric
    
    minimize ( y )
    subject to 
        (- Q0 + lam*H1 + y*H0) - NN == semidefinite(2*n+1) 
        NN(:) >= 0
cvx_end

ub = cvx_optval  ;  % from sdp we have obtained upper bound
utime = toc;

tic;
xl = [zeros(floor(n/2),1); ones(ceil(n/2),1)]; % just trivial feasible sol
% xl = bound.GW_lbound(2*[1 y'; y Y]+2,W(U,U),100);
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