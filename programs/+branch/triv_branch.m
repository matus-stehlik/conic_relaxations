function i = triv_branch(P,W)

n = size(W,1);
K = [P.pos,P.neg];  % known coordinates of x
if numel(K) == n, 
    i = -1;
    return;
end

for i = 1:n,                                            
   if isempty(find(P.pos == i,1)) && isempty(find(P.neg == i,1)),
       break;
   end
end


end