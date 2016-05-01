function xu = GW_ubound(Y, W, iter)

n = size(Y,1) ;
[U,p] = chol(Y);
if p>0,
    xu = bound.triv_bound(Y(2:end,2:end));
    return;
end


for i = 1:size(U,1),
    U(:,i) = U(:,i)/norm(U(:,i));
end

%a = 0;

ubound = inf;
for i = 1:iter,
    r = randn(n,1);
    r = r/norm(r);
    xu_new = U'*r;
    xu_new = xu_new(2:end);
    xu_new(xu_new >0) = 1;
    xu_new(xu_new <= 0) = -1;
    ubound_new = xu_new'*W*xu_new;
    
    %a = a*(i-1)/i + ubound_new/i;
    
    if ubound_new < ubound,
       ubound = ubound_new;
       xu = xu_new;
       
       %iTop = i;
    end
end

%fprintf(['iTop %d ,  average %4.3f ,    best found   %4.3f \n'],iTop,a,ubound)

end