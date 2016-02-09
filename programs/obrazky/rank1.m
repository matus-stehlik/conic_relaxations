% plot of symmetric rank 1 matrices under R^3
% the boundary is formed by matrices [x11,x12;x12,x22] 
% with 0 determinant i.e. x11*x22 - x12^2 = 0;  
% x11 = x12^2/x22 if x22~=0; and x11~=0 and x12=0 if x22=0;
n = 1500;
a = 1;
[x12,x22] = meshgrid(a/(2*n)*((-2*n):(2*n)));
x11 =  x12.^2./x22;

%x22(x22>2) = nan;
x11(abs(x11)>a*2) = nan;
    

mesh(x12,x22,x11, x11+x22);

axis equal
hold on;
mesh(x12,x11,x22,x11+x22);

hold off;
