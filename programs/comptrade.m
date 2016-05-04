% compasrison time-precision trade off

%% comp1
load('results1_final.mat');
Results = Results(11:end,[4:4:end,5:4:end,end]);

% %% comp3
% load('results2_60_80_100.mat')
% Results = Results2(21:end,[4:4:end,5:4:end,end]);


optVals = Results.optVal;

C = ...
 [ 0 0 1 % 1 BLUE
   0 1 0 % 2 GREEN (pale)
   1 0 0 % 3 RED
   0 1 1 % 4 CYAN
   1 0 1 % 5 MAGENTA (pale)
   1 1 0 % 6 YELLOW (pale)
   0 0 0 % 7 BLACK
   0 0.75 0.75 % 8 TURQUOISE
   0 0.5 0 % 9 GREEN (dark)
   0.75 0.75 0 % 10 YELLOW (dark)
   1 0.50 0.25 % 11 ORANGE
   0.75 0 0.75 % 12 MAGENTA (dark)
   0.7 0.7 0.7 % 13 GREY
   0.8 0.7 0.6 % 14 BROWN (pale)
   0.6 0.5 0.4 ]; % 15 BROWN (dark)

k = floor(width(Results)/2);

% upper bounds
rel_err = @(col)((col - optVals)./optVals);
rel_errors = varfun(rel_err,Results(:,1:k));
utimes = Results(:,(k+1):(end-1));

ub_names = cellstr(rel_errors.Properties.VariableNames);
ut_names = cellstr(utimes.Properties.VariableNames);
figure; hold on;
for i = 1:width(rel_errors)
    plot(rel_errors.(cell2mat(ub_names(i))),utimes.(cell2mat(ut_names(i))) ,'o','color', C(i,:), 'LineWidth',1.8 )
end
ax = gca;
ax.YScale = 'log';
ax.XScale = 'log';

names = ub_names;
for i = 1:length(names) 
    namei = cell2mat(names(i));
    names(i) = cellstr(regexprep(namei(5:(end-3)),'_',' '));
end

legend(names)
title('Comparison of relaxations on 10 rudy instances for n=80')
xlabel('precision (relative error)');
ylabel('running time');

% probably should also do 100 case and run those slow DNNs

