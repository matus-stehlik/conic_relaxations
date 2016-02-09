% dir1 = 'C:\Users\matusko\Documents\MATLAB\diplomka\problem_instances\biq\be' ;
% files = dir([dir1 '\*.mat']);
% m = length(files);
% optVals = zeros(m,1);
% for i = 1:m,
%     A = textread([dir1 '\' files(i).name]);
%     n = A(1,1);
%     A = A(2:end,:);
%     optVals(i) = sdp(A,n);
% end

dir1 = 'C:\Users\matusko\Documents\MATLAB\diplomka\problem_instances\mac\rudy' ;
files = dir([dir1 '\g05_60.*']);
m = length(files);
optVals = zeros(m,1);
for i = 1:1, %1:m 
    A = textread([dir1 '\' files(i).name]);
    n = A(1,1);    
    B = sparse(A(2:end,1),A(2:end,2),A(2:end,3),n,n);
    B = B + B' - 2*diag(diag(B)); %to make it symmetrical
    L = diag(B*ones(n,1)) - B;
    tic; [x,optVals(i), topIter] = bnb(-L*1/4); toc;
end


optVals = -optVals    
    