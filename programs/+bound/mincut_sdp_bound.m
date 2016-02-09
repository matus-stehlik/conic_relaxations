function [lb, ub, xf] = mincut_sdp_bound(W,P)

n = size(W,1);
ni = 1:n;

[mv, mw] = mincut2(W, P.pos, P.neg);

% let P.pos be in the 1 partition (not 0)
if all(mv(P.pos) == 0), mv = 1-mv; end 

P1 = P;
P1.neg = [];
for i = 1:length(P.pos)
    P1.pos(i) = P1.pos(i) - sum(mv(1:P1.pos(i))==0);
end

P0 = P;
P0.pos = [];
for i = 1:length(P.pos)
    P0.neg(i) = P0.neg(i) - sum(mv(1:P0.neg(i))==0);
end
% P0.pos = setdiff(P.pos, ni(mv==1));
% P0.neg = setdiff(P.neg, ni(mv==1));



L = diag(W*ones(n,1)) - W;

[lb1, ub1, xf1] = bound.sdp_bound(-1/4*L(mv==1,mv==1),P1);
[lb0, ub0, xf0] = bound.sdp_bound(-1/4*L(mv==0,mv==0),P0);


xf = zeros(n,1);
xf(mv==0) = xf0;

xfa = xf;
xf(mv==1) = xf1;

ub = -1/4*xf'*L*xf;


lb = lb1 + lb0 - mw ; %toto sa da zlepsit este o mincuty v ramci pos a neg



end
