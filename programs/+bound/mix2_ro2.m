function [lb, ub, xf, utime, ltime] = mix2_ro2(W,P)
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
cvx_solver sedumi

tic;
wd = diag(W(U,U));

% compute the largest eigenvalue 
lambda = eig(W(U,U)-diag(wd));
lam = lambda(end) + 1e-10; %correction to keep negative semidefinitnes in the comp. precision



cvx_begin
    cvx_quiet(true);
    variable x(n,1)
    maximize ( w0 + sum(wd) + lam*n + 2*wk'*x + x'*(W(U,U)-diag(wd) - lam*eye(n))*x )
    subject to 
         x.^2 <= 1        
cvx_end

ub = cvx_optval  ;  % from sdp we have obtained lower bound
utime = toc;

                            
% Rounding 1 - trivial boound.
tic;

alpha = 1 - x.^2;
beta = alpha;
ri = zeros(n,1);
ci = zeros(n,1);
vi = zeros(n,1);
count = 1;


[~,ordW] = sort(min(W(U,U)-diag(diag(W(U,U))) ));
for i = ordW(1:end-1), 
    if (beta(i)>0), 
        [col,indC] = sort(W(U(i),U([1:(i-1),(i+1):end])));
        for j = 1:length(col),
            if beta(indC(j))>0, 
                ri(count) = i;
                ci(count) = indC(j);
                if beta(indC(j))>beta(i),
                    vi(count) = beta(i)/2;
                    beta(indC(j)) = beta(indC(j)) - beta(i);
                    beta(i) = 0;
                    count = count + 1;
                    break;  
                else
                    vi(count) = beta(indC(j))/2;
                    beta(i) = beta(i) - beta(indC(j));
                    beta(indC(j)) = 0;
                    count = count +1;
                end
            end
        end
    end
end
ri = ri(vi~=0);
ci = ci(vi~=0);
vi = vi(vi~=0);
B = sparse([ri;ci],[ci;ri],[vi;vi],n,n);


% Y = [1, x'; x, (diag(alpha) - diag(alpha - ones(n,1)))*x*x'];
Y = [1, x'; x, x*x' + diag(alpha) - B];
xl = bound.GW_lbound(Y,W(U,U),1000);


xf = zeros(N,1);
xf(U) = xl;
xf(K) = xk;
lb1 = xl'*W(U,U)*xl + 2*wk'*xl+ w0 ;
lb2 = lb1 - 4*wk'*xl;
[lb,ind] = max([lb1,lb2]);
if ind==2, xf(U) = -xf(U); end
ltime = toc;
end