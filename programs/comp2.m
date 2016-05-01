%% Comparison 2

P = struct('pos', [1], 'neg', [], 'lbound', -inf , 'ubound', inf);

Results2 = table;

optVals60 = [536 532 529 538 527 533 531 535 530 533];
optVals80 = [929 941 934 923 932 926 929 929 925 923];
optVals100 = [1430 1425 1432 1424 1440 1436 1434 1431 1432 1430];
optVals = [optVals60 optVals80 optVals100]';

inst_len = 12;

cvx_clear;


dir1 = 'c:\Users\Miroslav\Desktop\diplomka\diplomka-master\programs\problem_instances\mac\rudy\' ;
files = [dir([dir1 'g05_60.*']); dir([dir1 'g05_80.*']); dir([dir1 'g05_100.*'])];
m = length(files);
for i = 1:m, %1:m 
    A = textread([dir1 '\' files(i).name]);
    n = A(1,1);    
    B = sparse(A(2:end,1),A(2:end,2),A(2:end,3),n,n);
    B = B + B' - 2*diag(diag(B)); %to make it symmetrical
    L = diag(B*ones(n,1)) - B;

    W = L/4;

    result = struct();
        
    result.inst = [files(i).name repmat(' ',1,inst_len-length(files(i).name))]; 
    result.Size = n;

    files(i).name

    [result.sdp2_lb, result.sdp2_ub, xf, result.sdp2_utime, result.sdp2_ltime] = bound.sdp2(W,P); 
    sprintf('sdp2 %d done',i)
    [result.sdp_lb, result.sdp_ub, xf, result.sdp_utime, result.sdp_ltime] = bound.sdp(W,P); 
    sprintf('sdp %d done',i)
%     [result.sdp_rlt_lb, result.sdp_rlt_ub, xf, result.sdp_rlt_utime, result.sdp_rlt_ltime] = bound.sdp_rlt(W,P);
%     sprintf('sdp_rlt %d done',i)
%     [result.sdp_tri_lb, result.sdp_tri_ub, xf, result.sdp_tri_utime, result.sdp_tri_ltime] = bound.sdp_triangle(W,P);
%     sprintf('sdp_tri %d done',i)
%     [result.socp1_lb, result.socp1_ub, xf, result.socp1_utime, result.socp1_ltime] = bound.socp1(W,P);
%     sprintf('socp1 %d done',i)
%     [result.socp2_lb, result.socp2_ub, xf, result.socp2_utime, result.socp2_ltime] = bound.socp2(W,P);
%     sprintf('socp2 %d done',i)
%     [result.socp3_lb, result.socp3_ub, xf, result.socp3_utime, result.socp3_ltime] = bound.socp3(W,P);
%     sprintf('socp3 %d done',i)
%     [result.socp4_lb, result.socp4_ub, xf, result.socp4_utime, result.socp4_ltime] = bound.socp4(W,P);
%     sprintf('socp4 %d done',i)
%     [result.dnn1_lb, result.dnn1_ub, xf, result.dnn1_utime, result.dnn1_ltime] = bound.dnn1(W,P);
%     sprintf('dnn1 %d done',i)
%     [result.dnn2_lb, result.dnn2_ub, xf, result.dnn2_utime, result.dnn2_ltime] = bound.dnn2(W,P);
%     sprintf('dnn2 %d done',i)
%     [result.dnn3p_lb, result.dnn3p_ub, xf, result.dnn3p_utime, result.dnn3p_ltime] = bound.dnn3p(W,P,10000);
%     sprintf('dnn3p %d done',i)
%     [result.dnn3d_lb, result.dnn3d_ub, xf, result.dnn3d_utime, result.dnn3d_ltime] = bound.dnn3d(W,P,10000);
%     sprintf('dnn3d %d done',i)
    [result.lp1_lb, result.lp1_ub, xf, result.lp1_utime, result.lp1_ltime] = bound.lp1(W,P);
    sprintf('lp1 %d done',i)
%     [result.lp2_lb, result.lp2_ub, xf, result.lp2_utime, result.lp2_ltime] = bound.lp2(W,P);
%     sprintf('lp2 %d done',i)
    [result.lp3_lb, result.lp3_ub, xf, result.lp3_utime, result.lp3_ltime] = bound.lp3(W,P);
    sprintf('lp3 %d done',i)
%     [result.lp4_lb, result.lp4_ub, xf, result.lp4_utime, result.lp4_ltime] = bound.lp4(W,P);
%     sprintf('lp4 %d done',i)
    [result.mix_lb, result.mix_ub, xf, result.mix_utime, result.mix_ltime] = bound.mixedSocpSdp(W,P);
    sprintf('mixr20 %d done',i)
    [result.mixr20_lb, result.mixr20_ub, xf, result.mixr20_utime, result.mixr20_ltime] = bound.mixedSocpSdpr(W,P,20);
    sprintf('mixr20 %d done',i)
    [result.mixri20_lb, result.mixri20_ub, xf, result.mixri20_utime, result.mixri20_ltime] = bound.mixedSocpSdpr2(W,P,20);
    sprintf('mixri20 %d done',i)
    [result.mixr3_lb, result.mixr3_ub, xf, result.mixr3_utime, result.mixr3_ltime] = bound.mixedSocpSdpr(W,P,3);
    sprintf('mixr3 %d done',i)
    [result.mixri3_lb, result.mixri3_ub, xf, result.mixri3_utime, result.mixri3_ltime] = bound.mixedSocpSdpr2(W,P,3);
    sprintf('mixri3 %d done',i)


Results2 = [Results2; struct2table(result)];
save('results2.mat', 'Results2');
    
    % started 10:50 PM
end

Results2.optVal = optVals;
Results2(:,[1,2, end,4:4:end])

%%


optVals60 = [536 532 529 538 527 533 531 535 530 533];
optVals80 = [929 941 934 923 932 926 929 929 925 923];
optVals100 = [1430 1425 1432 1424 1440 1436 1434 1431 1432 1430];
optVals = [optVals60 optVals80 optVals100]';



% upper bounds
utabl = Results2(:,[1,2, end,4:4:end]);
rel_err = @(col)((col - optVals)./optVals);
rel_errors = varfun(rel_err,utabl(:,4:end));

ub_names = cellstr(rel_errors.Properties.VariableNames);
figure; hold on;
for i = 1:width(rel_errors)
    plot(i,rel_errors.(cell2mat(ub_names(i))), 'ro')
end
ax = gca;
ax.YScale = 'log';




% times 
uttabl = Results2(:,[5:4:end]);
ut_names = cellstr(uttabl.Properties.VariableNames);
hold on;
for i = 1:width(uttabl)
    plot(i,uttabl.(cell2mat(ut_names(i))), 'b*')
end
ax = gca;
ax.YScale = 'log';


names = ub_names;
for i = 1:length(names) 
    namei = cell2mat(names(i));
    names(i) = cellstr(regexprep(namei(5:(end-3)),'_',' '));
end

set(gca,'xtick',1:length(names));
set(gca,'XTickLabel',names);
ax.XTickLabelRotation=90;
