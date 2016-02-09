% SDP relaxation of simple 0-1 programming: 
% min (something) 
% s.t.  xi^2 - xi = 0; for all i

% is problem:
% min ()
% s.t.  diag(X) = x 
%       X is symmetric 
%       X - xx^T is semidefinite

n = 10;
[x11,x22,x33,x12,x13,x23] = ndgrid(0:1/n:1,0:1/n:1,0:1/n:1,...
                                    -1:1/n:1,-1:1/n:1,-1:1/n:1);
% main subdeterminants must be nonegative
x11((x11-x11.^2).*(x22-x22.^2) - (x12-x11.*x22).^2 <0) = nan; 
x11(    (x11-x11.^2).*(x22-x22.^2).*(x33-x33.^2) +          ...
        2*(x12-x11.*x22).*(x23-x22.*x33).*(x13-x11.*x33) -  ...
        ((x13-x11.*x33).^2).*(x22-x22.^2) -                 ...
        (x11-x11.^2).*(x23-x22.*x33).^2 -                   ...
        (x33-x33.^2).*(x12-x11.*x22).^2 <0                  ) = nan;

    %sum(sum(sum(sum(sum(sum(isnan(x11)))))))  

%% 

x1 = x11(~isnan(x11));
x2 = x22(~isnan(x11));
x3 = x33(~isnan(x11));
figure(1);
plot3(x1,x2,x3,'r.');

k = convhull(x1,x2,x3);
%plot3(x1(k),x2(k),x3(k))

figure(2);
%dt =  DelaunayTri(x1(k),x2(k),x3(k));
trisurf(k,x1,x2,x3)
