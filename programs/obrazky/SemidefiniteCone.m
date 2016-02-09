% plot boundary of Semidefinite Cone S^2_+ as a cone under R^3
% the boundary is formed by matrices [x11,x12;x12,x22] 
% with 0 determinant i.e. x11*x22 - x12^2 = 0; and x11 >= 0; 
% x11 = x12^2/x22 for x22>0; and x11>=0 , x12=0 for x22=0;
n = 100;
a = 20;
% [x12,x22] = meshgrid(a/(2*n)*((-2*n):(2*n)), a/(n)*((1):(2*n)));
% x11 =  x12.^2./x22;
% % x11(x22+x11>2) = nan;
% % x22(x22+x11>2) = nan;

%x22(x22>2) = nan;
% x11(x11>a*2) = nan;

[x11,x22] = meshgrid(a/(n)*((1):(2*n)));
 x12 =  sqrt(x11.*x22);
 x11(x11+x22>a*2) = nan;
mesh(x12,x22,x11, x11+x22);

hold on;
mesh(-x12,x11,x22,x11+x22);

axis equal

hold off;
