fig1 = figure();
%shuffling_num = 10000;
axislabel_for_plot = {'pred vs unpred fam', 'Recency', 'Violation','Novelty','Reward value precue', 'Uncertainty precue split', 'RPE', 'Infoanticip split'};

% see script in Y:\Kaining\Infonew_analysis\index_infonew_code_testing
indices_for_plot = {'pred_vs_unpred_fam','recency_ind_match_pos','violation_ind', 'pred_nov_vs_fam',...
    'rewardvalueindex_precue', 'uncertaintyindex_precue_split', 'signedRPE_cue', 'RewInfoAnticipIndex_split'};
P_value_for_plot = {'P_pred_vs_unpred_fam_perm', 'P_recency_ind_match_pos','P_violation_ind_perm', 'P_pred_nov_vs_fam',...
    'rewardvalueindexP_precue','uncertaintyindexP_precue_split', 'signedRPE_cue_P', 'RewInfoAnticipIndexP_split'};

% axislabel_for_plot = {'Reward value precue', 'Uncertainty precue split', 'Simple RPE', 'Infoanticip split'};
% indices_for_plot = {'rewardvalueindex_precue', 'uncertaintyindex_precue_split','simpleRPErm_precueroc_cue','RewInfoAnticipIndex_split'};
% P_value_for_plot = {'rewardvalueindexP_precue', 'uncertaintyindexP_precue_split', 'simpleRPErm_precueroc_cue_P', 'RewInfoAnticipIndexP_split'};

    
% Pvalues_for_plot = {'Ppred_vs_unpred_fam','Precency_ind','Pviolation_ind'};
permutexy = {[4,5], [4,6],[4,7],[4,8], [1,5],[1,6],[1,7],[1,8],[2,5],[2,6],[2,7],[2,8] };

plotplacesetx = {11:40, 51:80,91:120,131:160,...
    11:40, 51:80,91:120,131:160,11:40, 51:80,91:120,131:160};
plotplacesety = {1:49, 1:49, 1:49, 1:49, 61:109, 61:109, 61:109, 61:109, 121:169,121:169,121:169, 121:169};
% barcolor = {'b','r','g'};
% Xedges = -0.5:0.05:0.5;
% Yedges = -0.5:0.05:0.5;

Include_neurons = find([Neuronlist_good(:).P_pred_nov_vs_fam]<StatisticalThreshold & [Neuronlist_good(:).pred_nov_vs_fam]>0 & ~isnan([Neuronlist_good(:).rewardvalueindex]));

%Include_neurons = find(~isnan([Neuronlist_good(:).rewardvalueindex]));

for xyw = 1:length(permutexy) %%% x-axis recency ind, y-axis violation_ind
    
%     if contains(indices_for_plot{permutexy{xyw}(1)}, {'rewardvalueindex','uncertaintyindex', 'RPE', 'RewInfoAnticipIndex'})
%         xaxis_ind =  [Neuronlist_good(:).(indices_for_plot{permutexy{xyw}(1)})]';
%     else
%         xaxis_ind =  [Neuronlist_good(:).(indices_for_plot{permutexy{xyw}(1)})]'-0.5;
%     end
%     
%     if contains(indices_for_plot{permutexy{xyw}(2)}, {'rewardvalueindex','uncertaintyindex', 'RPE', 'RewInfoAnticipIndex'})
%         yaxis_ind =  [Neuronlist_good(:).(indices_for_plot{permutexy{xyw}(2)})]';
%     else
%         yaxis_ind =  [Neuronlist_good(:).(indices_for_plot{permutexy{xyw}(2)})]'-0.5;
%     end
    xaxis_ind =  [Neuronlist_good(:).(indices_for_plot{permutexy{xyw}(1)})]';
    yaxis_ind =  [Neuronlist_good(:).(indices_for_plot{permutexy{xyw}(2)})]';
    
    xaxis_ind = xaxis_ind(Include_neurons);
    yaxis_ind = yaxis_ind(Include_neurons);
    
    NovelExcited_local=find([Neuronlist_good(Include_neurons).pred_nov_vs_fam]>0 & [Neuronlist_good(Include_neurons).P_pred_nov_vs_fam]<StatisticalThreshold)';
    NovelInhibited_local=find([Neuronlist_good(Include_neurons).pred_nov_vs_fam]<0 & [Neuronlist_good(Include_neurons).P_pred_nov_vs_fam]<StatisticalThreshold)';

    
    
    xaxis_label = axislabel_for_plot{permutexy{xyw}(1)};
    yaxis_label = axislabel_for_plot{permutexy{xyw}(2)};
    
    nsubplot(169,169, plotplacesety{xyw}, plotplacesetx{xyw}); set(gca,'ticklength',4*get(gca,'ticklength'))
    
    
    %% All neurons
    lim = [-1 1];
    nbin = 25;
    binedge = linspace(lim(1),lim(2),nbin+1);
    binmid = 0.5*binedge(1:(end-1)) + 0.5*binedge(2:end);
    
    %for xyx = 1 %1: non-novelty, 2: novelty excited, 3: novelty inhibited
    %eval(['currentindset = ' Indvariablenames{xyx} ';']);
    %histogram2(xaxis_ind(currentindset),yaxis_ind(currentindset), Xedges,  Yedges,'facecolor','flat');
    [n,xedge,yedge] = histcounts2(xaxis_ind,yaxis_ind,binedge,binedge);
    
    n = n'; % transpose so that x=columns and y=rows
    
    %n(isnan(n)) = 0;
    n = n ./ nansum(n(:));
    n = log(n);
    minimum_n = min(n(n>-inf));
    n(n==-inf)=minimum_n-0.5;
    colormap(flipud(gray));
    %figuren;
    %image(binmid,binmid,colormapify(n,[0 max(n(:))],'w',interpcolor('w','k',.5),'k','w'));
    imagesc(binmid,binmid,n);
    n_contour_levels = 6;
    contour(binmid,binmid,n,n_contour_levels,'linewidth',3);
    %end
    
    axis([-1 1 -1 1]); %%axis square;
    xlabel(xaxis_label)
    ylabel(yaxis_label)
    %zlabel('neuron count');
    line([-1 1],[0 0],'color',[0.3 0.3 0.3],'LineWidth',1);
    line([0 0],[-1 1],'color',[0.3 0.3 0.3],'LineWidth',1);
    %view(0,90);
    
    %%
    
    % now get the point which both x and y are not nan
    notnanlogic = ~isnan(xaxis_ind) & ~isnan(yaxis_ind);

    [slope,yint] = type_2_regression(xaxis_ind(notnanlogic), yaxis_ind(notnanlogic));
    pcs = pca([xaxis_ind(notnanlogic)-mean(xaxis_ind(notnanlogic)), yaxis_ind(notnanlogic) - mean(yaxis_ind(notnanlogic))]);
    line([-1 1], yint+slope.*[-1 1],'color','r','LineWidth',1);
    %[rho,rhop] = corr(x,y,'type','Spearman');
    
    % All datapoint value
    [rho,p] = corr(yaxis_ind , xaxis_ind, 'Type', 'Spearman');
    text(0.2,0.45, ['p all ' mat2str(p,4)])
    text(0.2,0.5, ['rho all ' mat2str(rho,4)])
    % positive novelty neuron
    try
        [rho_npos,p_npos] = corr(yaxis_ind(NovelExcited_local) , xaxis_ind(NovelExcited_local), 'Type', 'Spearman');
        text(0.2,0.35, ['p pos ' mat2str(p_npos,4)])
        text(0.2,0.4, ['rho pos ' mat2str(rho_npos,4)])
    catch
    end
    % negative novelty neuron
    try
    [rho_nneg,p_nneg] = corr(yaxis_ind(NovelInhibited_local) , xaxis_ind(NovelInhibited_local), 'Type', 'Spearman');
    text(0.2,0.25, ['p neg ' mat2str(p_nneg,4)])
    text(0.2,0.3, ['rho neg ' mat2str(rho_nneg,4)])
    catch
    end
    % non-related novelty neuron
    % [rho_nirl,p_nirl] = corr(yaxis_ind(NotNoveltySelective) , xaxis_ind(NotNoveltySelective), 'Type', 'Spearman');
    
    
    text(0.2,0.2, ['n = ' mat2str(sum(notnanlogic))])
%     text(0.2,0.15, ['p remain ' mat2str(p_nirl,4)])
%     text(0.2,0.2, ['rho remain ' mat2str(rho_nirl,4)])
    if strcmpi(indices_for_plot{permutexy{xyw}(1)}, 'pred_nov_vs_fam')
        xlim([-1 1]);
    end
    if strcmpi(indices_for_plot{permutexy{xyw}(2)}, 'pred_nov_vs_fam')
        ylim([-1 1]);
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make a table in which have rho number and indicated p value
fig2 = figure();
x_ind_name = {'Simple RPE', 'Reward value', 'Infoanticip'};
y_ind_name = {'Sensory', 'Recency'};
permutexy = {[1,7],[1,5],[1,8];
    [2,7],[2,5],[2,8]};
rhomatrix = zeros(size(permutexy)); % save rho number
pmatrix = zeros(size(permutexy)); % save p value
nmatrix = zeros(size(permutexy)); % save the number of the neuron
for ii = 1:size(permutexy,1)
    for jj = 1:size(permutexy,2)
        
%         if contains(indices_for_plot{permutexy{ii,jj}(1)}, {'rewardvalueindex','uncertaintyindex', 'RPE', 'RewInfoAnticipIndex'})
%             xaxis_ind =  [Neuronlist_good(:).(indices_for_plot{permutexy{ii,jj}(1)})]';
%         else
%             xaxis_ind =  [Neuronlist_good(:).(indices_for_plot{permutexy{ii,jj}(1)})]'-0.5;
%         end
%         
%         if contains(indices_for_plot{permutexy{ii,jj}(2)}, {'rewardvalueindex','uncertaintyindex', 'RPE', 'RewInfoAnticipIndex'})
%             yaxis_ind =  [Neuronlist_good(:).(indices_for_plot{permutexy{ii,jj}(2)})]';
%         else
%             yaxis_ind =  [Neuronlist_good(:).(indices_for_plot{permutexy{ii,jj}(2)})]'-0.5;
%         end
        xaxis_ind =  [Neuronlist_good(:).(indices_for_plot{permutexy{ii,jj}(1)})]';
        yaxis_ind =  [Neuronlist_good(:).(indices_for_plot{permutexy{ii,jj}(2)})]';
        
        xaxis_ind = xaxis_ind(Include_neurons);
        yaxis_ind = yaxis_ind(Include_neurons);
        
        notnanlogic = ~isnan(xaxis_ind) & ~isnan(yaxis_ind);
        [rho,p] = corr(yaxis_ind(notnanlogic) , xaxis_ind(notnanlogic), 'Type', 'Spearman');
        rhomatrix(ii,jj)= rho;
        pmatrix(ii,jj) = p;
        nmatrix(ii,jj) = sum(notnanlogic);
        
    end
end

heatmap(x_ind_name, y_ind_name, rhomatrix, 'Colormap', gray);
%text(1,1,'***')
%text(perc(segmentInd) + adj/2, barposition(segmentInd) + w/3,'***','Color',color_star);

print(gcf,'-dpdf', '-painters',fullfile(plotpath,['correlation_table_reward' '.pdf']));


%% a mega figure summarizing the info task's finding.
figure;
% show that 3 reward related indices's measure is above chance
nsubplot(169,169, plotplacesety{1}, plotplacesetx{1}); set(gca,'ticklength',4*get(gca,'ticklength'))
Pvaluenames = {'rewardvalueindexP_precue','signedRPE_cue_P', 'RewInfoAnticipIndexP_split'};
Pvaluelabel = {'Reward value', 'Reward surprise', 'Infoanti'};

sig_percent_vals = zeros(size(Pvaluenames));
p_binomial_vals = zeros(size(Pvaluenames));
error_bar_vals = zeros(size(Pvaluenames));
for xyw = 1:numel(Pvaluenames)
    Pvalues = [Neuronlist_good(:).(Pvaluenames{xyw})];
    Pvalues = Pvalues(~isnan(Pvalues));
    sig_num = sum(Pvalues<StatisticalThreshold);
    total_num = sum(~isnan(Pvalues));
    sig_percent_vals(xyw) = sig_num/total_num;
    Sided = 'Greater';
    p_binomial_vals(xyw) = myBinomTest(sig_num,total_num,StatisticalThreshold,Sided);
    error_bar_vals(xyw) = std(double(Pvalues<StatisticalThreshold))/sqrt(total_num);
end

bar([1:numel(Pvaluenames)], sig_percent_vals);
errorbar([1:numel(Pvaluenames)], sig_percent_vals, error_bar_vals);
plot(xlim, StatisticalThreshold*[1,1], 'b--');

for xyw = 1:numel(Pvaluenames)
    text(xyw,0.05, sprintf('p=%.3f',p_binomial_vals(xyw)));
end
set(gca, 'xtick', [1:numel(Pvaluenames)], 'xticklabel', Pvaluelabel);
ylabel('Percentage of significant neurons');
% for ii = 1:numel(segmentNames)
%     % write the actual #
%     % sig of binomial test
%     if P_corr.(currenttype)(segmentInd,condInd) <= pars.tableColor_low
%         text(perc(segmentInd) + adj/2, barposition(segmentInd) + w/3,'***','Color',color_star);
%     elseif P_corr.(currenttype)(segmentInd,condInd) <= pars.tableColor_mid
%         text(perc(segmentInd) + adj/2, barposition(segmentInd) + w/3,'**','Color',color_star);
%     elseif P_corr.(currenttype)(segmentInd,condInd) <= pars.tableColor_high
%         text(perc(segmentInd) + adj/2, barposition(segmentInd) + w/3,'*','Color',color_star);
%     end
%     %
% end

% a bar plot summarizing the correlation of correlation of Novelty index with the indices in NFL task and
% info tasks.

corr_permutexy = {[4,1], [4,2], [4,5], [4,7],[4,8], [1,5], [1,7],[1,8], [2,5], [2,7],[2,8]};
corr_label = {'Sensory surprise', 'Recency', 'Reward value', 'Reward surprise', 'Infoanti'...
    'Reward value', 'Reward surprise', 'Infoanti', 'Reward value', 'Reward surprise', 'Infoanti'};

corrvalue = zeros(size(corr_permutexy));
corrp = zeros(size(corr_permutexy));
corrstd = zeros(size(corr_permutexy));
corr_n = zeros(size(corr_permutexy));
for xyw = 1: numel(corr_permutexy)
%     if contains(indices_for_plot{corr_permutexy{xyw}(1)}, {'rewardvalueindex','uncertaintyindex', 'RPE', 'RewInfoAnticipIndex'})
%         xaxis_ind =  [Neuronlist_good(:).(indices_for_plot{corr_permutexy{xyw}(1)})]';
%     else
%         xaxis_ind =  [Neuronlist_good(:).(indices_for_plot{corr_permutexy{xyw}(1)})]'-0.5;
%     end
%     
%     if contains(indices_for_plot{corr_permutexy{xyw}(2)}, {'rewardvalueindex','uncertaintyindex', 'RPE', 'RewInfoAnticipIndex'})
%         yaxis_ind =  [Neuronlist_good(:).(indices_for_plot{corr_permutexy{xyw}(2)})]';
%     else
%         yaxis_ind =  [Neuronlist_good(:).(indices_for_plot{corr_permutexy{xyw}(2)})]'-0.5;
%     end
    
    xaxis_ind =  [Neuronlist_good(:).(indices_for_plot{corr_permutexy{xyw}(1)})]';
    yaxis_ind =  [Neuronlist_good(:).(indices_for_plot{corr_permutexy{xyw}(2)})]';
    
    xaxis_ind = xaxis_ind(Include_neurons);
    yaxis_ind = yaxis_ind(Include_neurons);
    
    notnanlogic = ~isnan(xaxis_ind) & ~isnan(yaxis_ind);
    [rho,p] = corr(xaxis_ind(notnanlogic) , yaxis_ind(notnanlogic), 'Type', 'Spearman');
    corrvalue(xyw) = rho;
    corrp(xyw) = p;
    corr_n(xyw) = numel(xaxis_ind);
    
    % bootstrapping to get std interval
    xaxis_ind = xaxis_ind(notnanlogic);
    yaxis_ind = yaxis_ind(notnanlogic);
    
    corr_shuffled = zeros(shuffling_num,1);
    for ii = 1:shuffling_num
        shuffling_ind = randi(numel(xaxis_ind), size(xaxis_ind));
        corr_shuffled(ii) = corr(xaxis_ind(shuffling_ind) , yaxis_ind(shuffling_ind), 'Type', 'Spearman');
    end
    corrstd(xyw) = std(corr_shuffled);
    
    
end

%plots
nsubplot(169,169, plotplacesety{2}, plotplacesetx{2}); set(gca,'ticklength',4*get(gca,'ticklength'))
bar([1:5],corrvalue(1:5));
errorbar([1:5],corrvalue(1:5), corrstd(1:5));
ylabel('correlation with novelty index');
ylim([-0.05, 0.15]);
set(gca, 'xtick', [1:5], 'xticklabel', corr_label(1:5));

nsubplot(169,169, plotplacesety{3}, plotplacesetx{3}); set(gca,'ticklength',4*get(gca,'ticklength'))
bar([1:3],corrvalue(5+(1:3)));
errorbar([1:3],corrvalue(5+(1:3)), corrstd(5+(1:3)));
ylabel('correlation with sensory surprise index');
ylim([-0.05, 0.15]);
set(gca, 'xtick', [1:3], 'xticklabel', corr_label(5+(1:3)));

nsubplot(169,169, plotplacesety{4}, plotplacesetx{4}); set(gca,'ticklength',4*get(gca,'ticklength'))
bar([1:3],corrvalue(8+(1:3)));
errorbar([1:3],corrvalue(8+(1:3)), corrstd(8+(1:3)));
ylabel('correlation with recency index');
ylim([-0.05, 0.15]);
set(gca, 'xtick', [1:3], 'xticklabel', corr_label(8+(1:3)));

% bar plot of correlation of recency and sensory surprise indices


set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
print(gcf,'-dpdf', '-painters',fullfile(plotpath,['NLF_reward_task_corr_comp_barplot' '.pdf']));



% save the plot
figure(fig1);
set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
print(gcf,'-dpdf', '-painters',fullfile(plotpath,['Indices_heatmap_withreward_all_session' '.pdf']));
% close gcf;


