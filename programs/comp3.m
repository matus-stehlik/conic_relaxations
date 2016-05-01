%% Comparison 3

P = struct('pos', [1], 'neg', [], 'lbound', -inf , 'ubound', inf);

Results3 = table;
cvx_clear;

m = 4;

for n = 100:100:600
    for i = 1:m, %1:m 
         B = round(rand(n)/1.5).*round(rand(n)*25);
         B = B' + B;
         B = B - diag(diag(B))  ;
         L = diag(B*ones(n,1)) -  B;


        W = L/4;

        result = struct();
        result.Size = n;
        

    %     [result.sdp2_lb, result.sdp2_ub, xf, result.sdp2_utime, result.sdp2_ltime] = bound.sdp2(W,P); 
    %     sprintf('sdp2 %d.%d done',[n,i])
        
        bound.sdp(ones(100),P); % warm up so that processor is already preparer
        sprintf('mock_sdp %d.%d done',[n,i])
        
        [result.sdp_lb, result.sdp_ub, xf, result.sdp_utime, result.sdp_ltime] = bound.sdp(W,P); 
        sprintf('sdp %d.%d done',[n,i])
        [result.lp1_lb, result.lp1_ub, xf, result.lp1_utime, result.lp1_ltime] = bound.lp1(W,P);
        sprintf('lp1 %d.%d done',[n,i])
        [result.lp3_lb, result.lp3_ub, xf, result.lp3_utime, result.lp3_ltime] = bound.lp3(W,P);
        sprintf('lp3 %d.%d done',[n,i])
        [result.mix1_lb, result.mix1_ub, xf, result.mix1_utime, result.mix1_ltime] = bound.mixedSocpSdp(W,P);
        sprintf('mix1 %d.%d done',[n,i])
        [result.mix2_lb, result.mix2_ub, xf, result.mix2_utime, result.mix2_ltime] = bound.mixedSocpSdp2(W,P);
        sprintf('mix2 %d.%d done',[n,i])
        [result.mixr50_lb, result.mixr50_ub, xf, result.mixr50_utime, result.mixr50_ltime] = bound.mixedSocpSdpr(W,P,50);
        sprintf('mixr50 %d.%d done',[n,i])
        [result.mixr20_lb, result.mixr20_ub, xf, result.mixr20_utime, result.mixr20_ltime] = bound.mixedSocpSdpr(W,P,20);
        sprintf('mixr20 %d.%d done',[n,i])
        [result.mixr10_lb, result.mixr10_ub, xf, result.mixr10_utime, result.mixr10_ltime] = bound.mixedSocpSdpr(W,P,10);
        sprintf('mixr10 %d.%d done',[n,i])
        [result.mixr5_lb, result.mixr5_ub, xf, result.mixr5_utime, result.mixr5_ltime] = bound.mixedSocpSdpr(W,P,5);
        sprintf('mixr5 %d.%d done',[n,i])



    Results3 = [Results3; struct2table(result)];
    save('results3.mat', 'Results3');

    end
end

%Results3.optVal = optVals;
Results3(:,[1,2, end,4:4:end])

%%

% upper bounds
utabl = Results3(:,[1,3:4:end]);
best_vals = min(utabl{:,2:end}')';
rel_err = @(col)((col - best_vals)./best_vals);
rel_errors = varfun(rel_err,utabl(:,2:end));

ub_names = cellstr(rel_errors.Properties.VariableNames);
figure; hold on;
for i = 1:width(rel_errors)
    plot(i,rel_errors.(cell2mat(ub_names(i))), 'ro')
end
ax = gca;
ax.YScale = 'log';




% times 
uttabl = Results3(:,[4:4:end]);
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

ro = plot(-1,-1,'ro');
bs = plot(-1,-1,'b*');
legend([ro,bs],'relative lower bound', 'runing time', 'Location','north')



% probably should also do 100 case and run those slow DNNs

