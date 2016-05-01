function xu = triv_bound(X)
%TRIV_UBOUND Summary of this function goes here
%   Detailed explanation goes here

xu = sign(X(:,1));    % some feasible solution
xu(xu == 0) = 1;

end

