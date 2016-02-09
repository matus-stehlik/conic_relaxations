function [x,optVal] = bnb_alloc(W)
% input: W symmetric n × n matrix
% output: x in {-1, 1}^n optimal solution, optVal optimal value

% initialize:
n = size(W,1);
Q(n^2) = struct('pos', [], 'neg', [], 'lbound', [] ); %allocate space

%zlb = lower_bound(W,[1],[]);
xbk = [ones(floor(n/2),1);zeros(ceil(n/2),1)];          % initial cut vector
zbk = xbk'*W*xbk;                                       % best known (bk) solution value
Q(1) = struct('pos', [1], 'neg', [], 'lbound', -inf );  % problem list
Q(2) = struct('pos', [], 'neg', [1], 'lbound', -inf );
NoProblems = 2; % number of problems in Q
FirstProblem = 1;
LastProblem = 2;
Iter = 0;


while NoProblems,
Iter = Iter + 1;    
% remove problem P from Q having zlb minimal
% fprintf('aktualne Q \n');
% for i = 1:NoProblems
%     Q(i)
% end
% fprintf('\n');

P = Q(1) ;
FirstProblem = 2;    %Q(1) = [];
NoProblems = NoProblems -1;
[P.lbound,ubound, xf] = sdp_bound(W,P);
if ubound < zbk,
    zbk = ubound;
    xbk = xf;
end

% cut problems with lbound > ubound
if NoProblems > 0,
    for i = NoProblems:-1:1, 
        if Q(i).lbound > zbk
            Q(i) = [];
            NoProblems = NoProblems-1;
            LastProblem = LastProblem-1;
        else break;
        end
    end  
end


% we cannot branch if whole x is given
if length(P.pos) + length(P.neg) == n,                  
    continue;
end

% determine xi to branch for P              
                                %THIS should be done in more clever way!!!
for i = 1:n,                                            
   if isempty(find(P.pos == i,1)) && isempty(find(P.neg == i,1)),
       break;
   end
end

% branch and add new problems to Q
P1 = P;                                                 
P1.neg(end+1) = i;
P.pos(end+1) = i;
Q(LastProblem+1) = P1;
LastProblem = LastProblem+1;
Q(LastProblem+1) = P;
LastProblem = LastProblem+1;    
NoProblems = NoProblems + 2;


% zotriedime list problémov
Q(1:NoProblems) = Q(FirstProblem:LastProblem);
FirstProblem = 1;
LastProblem = NoProblems;
[temp, ind] = sort([Q(FirstProblem:LastProblem).lbound]);
Q(FirstProblem:LastProblem) = Q(ind);


end %while

% return solution
x = xbk;
optVal = zbk;
Iter;

return;

end