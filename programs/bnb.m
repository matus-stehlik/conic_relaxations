function [x,optVal, topIter] = bnb(W)
% input: W symmetric n × n matrix
% output: x in {-1, 1}^n optimal solution, optVal optimal value

% initialize:
n = size(W,1);
%zlb = lower_bound(W,[1],[]);
xbk = [ones(floor(n/2),1);zeros(ceil(n/2),1)];          % initial cut vector
zbk = xbk'*W*xbk;                                       % best known (bk) solution value
Q(1) = struct('pos', [1], 'neg', [2], 'lbound', -inf , 'ubound', inf);  % problem list
Q(2) = struct('pos', [1 2], 'neg', [], 'lbound', -inf, 'ubound', inf);    
NoProblems = 2; % number of problems in Q
Iter = 0;
maxIter = 1000;
ConsidProblems = ones(maxIter,1);
excluded = 0;

figure(1);
hold on;

while ~isempty(Q) && Iter<maxIter,
Iter = Iter + 1;   
if mod(Iter,20)== 0, 
    Best_Known_so_far = zbk
    NoProblems
end
% remove problem P from Q having zlb minimal
P = Q(1)
Q(1) = [];
NoProblems = NoProblems -1;

% we cannot branch if whole x is given
if length(P.pos) + length(P.neg) == n,                  
    continue;
end

% determine xi to branch for P              
                                                         
%i = branch.triv_branch(P,W);
i = branch.w_branch(P,W);
if i ~= -1, %we can branch

    % branch and add new problems to Q
    P1 = P;                                                 
    P1.neg(end+1) = i;
    P.pos(end+1) = i;

    % it may be possible to optimize, since P and P1 differ only in 1
    % fixed variable 
    [P.lbound,P.ubound, xf] = bound.sdp_bound(W,P);
    if P.ubound < zbk,
        zbk = P.ubound
        xbk = xf;
        topIter = Iter;
    end

    [P1.lbound,P1.ubound, xf] = bound.sdp_bound(W,P1);
    if P1.ubound < zbk,
        zbk = P1.ubound
        xbk = xf;
        topIter = Iter;
    end

    % cut problems with lbound > ubound
    if NoProblems > 0,
        
        for i = NoProblems:-1:1, 
            if Q(i).lbound > zbk
                excluded = excluded + 2^(1 - numel(Q(i).pos) - numel(Q(i).neg));
                Q(i) = [];
                NoProblems = NoProblems-1;
            else break;
            end
        end 
    end
    if P1.lbound < zbk,
        Q(end+1) = P1;
        NoProblems = NoProblems + 1;
    else
        excluded = excluded + 2^(1 - length(P1.pos) - length(P1.neg));
    end
    if P.lbound < zbk,
        Q(end+1) = P;
        NoProblems = NoProblems + 1;
    else
        excluded = excluded + 2^(1 - length(P.pos) - length(P.neg));
    end

    % zotriedime list problémov
    [temp, ind] = sort([Q.lbound]);
    Q = Q(ind);
    
    ConsidProblems(Iter) = ConsidProblems(Iter) - excluded;
    
    plot(Iter,P.ubound, 'g.');
    plot(Iter,P.lbound, 'b.');   

    plot(Iter,P1.ubound, 'g.');
    plot(Iter,P1.lbound, 'b.');
    
    plot(Iter,zbk,'r*');
    
    
end
end %while

hold off
figure;
plot(1:Iter,ConsidProblems(1:Iter));

% return solution
x = xbk;
optVal = zbk;
topIter
Iter

return;

end