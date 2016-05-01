function xl = GW_lbound(Y, W, iter)

n = size(Y,1) ;
[U,p] = chol(Y);
if p>0,
    xl = bound.triv_bound(Y(2:end,2:end));
    return;
end


for i = 1:size(U,1),
    U(:,i) = U(:,i)/norm(U(:,i));
end

%a = 0;

lbound = -inf;
for i = 1:iter,
    r = randn(n,1);
    r = r/norm(r);
    xl_new = U'*r;
    xl_new = xl_new(2:end);
    xl_new(xl_new >0) = 1;
    xl_new(xl_new <= 0) = -1;
    lbound_new = xl_new'*W*xl_new;
    
    %a = a*(i-1)/i + ubound_new/i;
    
    if lbound_new > lbound,
       lbound = lbound_new;
       xl = xl_new;
       
       %iTop = i;
    end
end

%fprintf(['iTop %d ,  average %4.3f ,    best found   %4.3f \n'],iTop,a,ubound)

end