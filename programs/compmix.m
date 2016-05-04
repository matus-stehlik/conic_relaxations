%% comparison mix
load('results3_final.mat')
Results = Results3; 
%%

C = {'k','b','r','g','m',[.5 .6 .7],[.8 .2 .6]} % Cell array of colros.

% upper bounds
utabl = Results(:,[1,3,15:4:end]);
best_vals = min(utabl{:,2:end}')';
rel_err = @(col)((col - best_vals)./best_vals);
rel_errors =  varfun(rel_err,utabl(:,2:end));
rel_errors.Size = Results.Size;


ub_names = cellstr(rel_errors.Properties.VariableNames);
figure; hold on;

        %%just for the legend purposes 
        for i = 1:(width(rel_errors)-1) 
             plot(-1.-1,  's', 'color', cell2mat(C(i)) , 'MarkerFaceColor',cell2mat(C(i)));
        end
        
for i = 1:(width(rel_errors)-1)
     plot(rel_errors.Size, rel_errors.(cell2mat(ub_names(i))),  'o', 'color', cell2mat(C(i)),'LineWidth',1 )
     mean_vals = mean(reshape(rel_errors.((cell2mat(ub_names(i)))),4,6));
     plot(Results.Size(1:4:end), mean_vals,  'o-', 'color', cell2mat(C(i)),'LineWidth',0.8 )
end
ax = gca;
ax.YScale = 'log';




% times 
uttabl = Results(:,[4,16:4:end]);
ut_names = cellstr(uttabl.Properties.VariableNames);
hold on;
for i = 1:width(uttabl)
    plot(Results.Size, uttabl.(cell2mat(ut_names(i))),  '*', 'color', cell2mat(C(i)), 'LineWidth',1 )
    mean_vals = mean(reshape(uttabl.(cell2mat(ut_names(i))),4,6));
    plot(Results.Size(1:4:end), mean_vals,  '*-', 'color', cell2mat(C(i)),'LineWidth',0.8 )
end
ax = gca;
ax.YScale = 'log';


names = ub_names;
for i = 1:length(names) 
    namei = cell2mat(names(i));
    names(i) = cellstr(regexprep(namei(5:(end-3)),'_',' '));
end

set(gca,'xtick',100:100:600);
%set(gca,'XTickLabel',names);

legend(names(1:end-1), 'Location', 'northwest')

% ro = plot(-1,-1,'ro');
% bs = plot(-1,-1,'b*');
% legend([ro,bs],'relative lower bound', 'runing time', 'Location','north')
title('Mixed SOCP-SDP relaxations, effect of increasing n')
xlabel('n');
ylabel('relative error   |    running time             ');

% probably should also do 100 case and run those slow DNNs

