function xu = triv_ubound(X)
%TRIV_UBOUND Summary of this function goes here
%   Detailed explanation goes here

xu = X(:,1);    % some feasible solution
xu(xu>0) = 1;
xu(xu<=0) = -1;

end

