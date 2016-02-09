% plot boundary of Second Order Cone x1 >= ||(x2,x3)||.
n = 200;
[y,z] = meshgrid(1/(2*n)*((-2*n):(2*n)));
x = sqrt(y.^2 + z.^2);
x(x>1) = NaN; % cut it at height 1
mesh(y,z,x) % plot
axis equal;