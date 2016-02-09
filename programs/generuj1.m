function L = generuj1(n,max_weight,density)
% generate Laplacian of the size n integer max cut problem with 
% uniformly random weights from 0 to max_weight and given density

W = round(triu((sprand(n,n,density)),1)*max_weight);
W = W' + W;
L = diag(W*ones(n,1)) - W;


return;

end