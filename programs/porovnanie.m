% porovnanie sdp_bound pre problem W,P a mincut_sdp_bound

% n = 100;
% W = magic(n);
% W = W'+W;


dir1 = 'C:\Users\matusko\Documents\MATLAB\diplomka\problem_instances\mac\rudy' ;
files = dir([dir1 '\p*']);
m = length(files);
optVals = zeros(m,1);
P = struct('pos',1:10, 'neg', 11:20);
for i = 1:m, %1:m 
    A = textread([dir1 '\' files(i).name]);
    n = A(1,1);    
    B = sparse(A(2:end,1),A(2:end,2),A(2:end,3),n,n);
    W = B + B' - 2*diag(diag(B)); %to make it symmetrical
    L = diag(W*ones(n,1)) - W;



    tic; [lb0, ub0, xf0] = bound.mincut_sdp_bound(W,P); toc;
    tic; [lb1, ub1, xf1] = bound.sdp_bound(-1/4*L,P); toc;

    rozdiely = [lb1-lb0, ub0-ub1]


    [lb0,ub0; lb1, ub1]
    
end
    
