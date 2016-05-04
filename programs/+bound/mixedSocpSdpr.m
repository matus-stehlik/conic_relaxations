function [lb, ub, xf, utime, ltime] = mixedSocpSdpr(W,P,r)
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
cvx_solver sedumi

% if all the variables are set, we dont need any optimization
if numel(K) == N,
    utime = -1;
    ltime = -1;
    lb = w0;
    ub = w0;
    xf = xk;
    return;
end


% Find a right C-block diagonal matrix based on the article
% by Burer, Kim and Kojima, Faster, but Weaker, Relaxations for
% Quadratically Constrained Quadratic Programs

% r is given number of blocks, we will choose the blocks uniformly
% 1,...,floor(n/r) | floor(n/r)+1,...,floor(2n/r) | ... |
% floor((r-1)n/r)+1,...,n


poc2 = mod(n,r);
poc1 = r-poc2;
m = floor(n/r);

C1 = zeros(m,poc1);
C2 = zeros(m+1,poc2);
count1 = 0;
count2 = 0;
for i = 1:r
    Ci = (floor((i-1)*n/r)+1):floor(i*n/r);
    if length(Ci) == m
        count1 = count1 + 1;
        C1(:,count1) = Ci;
    else 
        count2 = count2 + 1;
        C2(:,count2) = Ci;
    end
end

A = W(U,U);

tic;
for i = 1:r,
    A((floor((i-1)*n/r)+1):floor(i*n/r),(floor((i-1)*n/r)+1):floor(i*n/r)) = 0;    
end
lambda = eig(A);
lam = lambda(end)+ 1e-10; %correction to keep B negative semidefinite
B = A - lam*eye(n)-eye(n); 

D = sparse(W(U,U)-B); % the block diagonal matrix by the second shift of A


% min x^TWx with constraints on x given in P is equivalent to 
% min xu^TW(U,U)xu + 2wk^Tx + w0, s.t. xu is from {-1,1}^(n-1)
% which can be simplified to min trace(MY), s.t. diag(Y)=1 and Y is PSD
% where M = [w0,wk';wk,W(U,U)]

cvx_begin
    cvx_quiet(true);
    variable x(n,1)
    variable X(n,n) symmetric 
   variable S(r,1)
    maximize ( w0 + 2*wk'*x + sum(S) + x'*(B)*x ) 
    subject to 
         x.^2 <= 1
         diag(X) == 1
         X(D==0) == 0
         for i = 1:poc1
            [1 x(C1(:,i))' ;x(C1(:,i)) X(C1(:,i),C1(:,i))] == semidefinite(m+1) 
             S(i) == trace(X(C1(:,i),C1(:,i))*D(C1(:,i),C1(:,i)))
         end
         for i = 1:poc2
            [1 x(C2(:,i))'; x(C2(:,i)) X(C2(:,i),C2(:,i))] == semidefinite(m+2)
             S(i+poc1) == trace(X(C2(:,i),C2(:,i))*D(C2(:,i),C2(:,i)))
         end
cvx_end

ub = cvx_optval  ;  % from sdp we have obtained lower bound
utime = toc;
             
tic;
% !!! this should be done in more clever way
xl = bound.triv_bound(x);

xf = zeros(N,1);
xf(U) = xl;
xf(K) = xk;
lb1 = xl'*W(U,U)*xl + 2*wk'*xl+ w0 ;
lb2 = lb1 - 4*wk'*xl;
[lb,ind] = max([lb1,lb2]);
if ind==2, xf(U) = -xf(U); end
ltime = toc;

end