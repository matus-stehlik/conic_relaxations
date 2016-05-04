%% comparison of Branch and Bound perfomance of Mix2 and SDP

n = 40;
m = 3;
density = 0.5;
%ResultsBnb = table;
load('results_bnb_fin.mat');


for i = 1:m 
    
    W = (1/4) * generate_MC(n,density,0,25);
    result = struct();
    
    result.size = n;
    result.density = density;
    
%     tic;
%     [x,result.opt_sdp_slow, result.topIter_sdp_slow, result.Iter_sdp_slow] = bnb(W, @bound.sdp_slow_bnb);
%     result.time_sdp_slow = toc;
%     sprintf('sdp_slow %d.%d done',[n,i])
%     
    tic;
    [x,result.opt_sdp, result.topIter_sdp, result.Iter_sdp] = bnb(W, @bound.sdp_bnb);
    result.time_sdp = toc;
    sprintf('sdp %d.%d done',[n,i])
    
    tic;
    [x,result.opt_mix, result.topIter_mix, result.Iter_mix] = bnb(W, @bound.mix2_bnb);
    result.time_mix = toc;
    sprintf('mix %d.%d done',[n,i])
    

    
    ResultsBnb = [ResultsBnb; struct2table(result)];
    save('results_bnb_fin.mat', 'ResultsBnb');
    
end
    %%
    
    % considered problems
    legend('sdp', 'mix')
    title(sprintf('Branch & Bound for max cut with n = %d, density = %1.2f', ...
        [n,density]));
    xlabel('iteration')
    ylabel('portion of problems being considered');
    
    
    
    % SDP
    title(sprintf('Branch & Bound for max cut with SDP , n = %d, density = %1.2f', ...
            [n,density]));
    xlabel('iteration')
    ylabel('lower bound      |     upper bound');
    ti = result.Iter_sdp %Iter;
    ov = result.opt_sdp %optVal;
    
    hold on;
    p1 = plot(ti,ov,'b.');
    p2 = plot(ti,ov,'g.');
    p3 = plot(ti,ov,'r*');
    legend([p2,p3,p1], 'upper bound', 'best known solution', 'lower bound')
    hold off;
    
    
    
    %MIX
    title(sprintf('Branch & Bound for max cut with Mix , n = %d, density = %1.2f', ...
            [n,density]));
    xlabel('iteration')
    ylabel('lower bound      |     upper bound');
    
        ti = result.Iter_mix %Iter;
    ov = result.opt_mix %optVal;
    
    hold on;
    p1 = plot(ti,ov,'b.');
    p2 = plot(ti,ov,'g.');
    p3 = plot(ti,ov,'r*');
    legend([p2,p3,p1], 'upper bound', 'best known solution', 'lower bound')
    hold off;
    

    