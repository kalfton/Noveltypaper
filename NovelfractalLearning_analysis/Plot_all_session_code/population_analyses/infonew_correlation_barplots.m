fig1 = figure();
%shuffling_num = 10000;
axislabel_for_plot = {'Surprise', 'Recency', 'Violation','Novelty','Reward value', 'Infoanticip'};

% see script in Y:\Kaining\Infonew_analysis\index_infonew_code_testing
indices_for_plot = {'pred_vs_unpred_fam','recency_ind_match_pos','violation_ind', 'pred_nov_vs_fam',...
    'rewardvalueindex_precue', 'RewInfoAnticipIndex_split'};
P_value_for_plot = {'P_pred_vs_unpred_fam_perm', 'P_recency_ind_match_pos','P_violation_ind_perm', 'P_pred_nov_vs_fam',...
    'rewardvalueindexP_precue','RewInfoAnticipIndexP_split'};

plotplacesetx = {11:40, 51:80,91:120,131:160,...
    11:40, 51:80,91:120,131:160,11:40, 51:80,91:120,131:160};
plotplacesety = {1:49, 1:49, 1:49, 1:49, 61:109, 61:109, 61:109, 61:109, 121:169,121:169,121:169, 121:169};


Include_type = 'Noveltyexcited'; %'Noveltyexcited' or 'Noveltyinhibited';% 
if strcmpi(Include_type, 'Noveltyexcited')
    Include_neurons = find([Neuronlist_good(:).P_pred_nov_vs_fam]<StatisticalThreshold & [Neuronlist_good(:).pred_nov_vs_fam]>0 & ~isnan([Neuronlist_good(:).rewardvalueindex_precue]));
elseif strcmpi(Include_type, 'Noveltyinhibited')
    Include_neurons = find([Neuronlist_good(:).P_pred_nov_vs_fam]<StatisticalThreshold & [Neuronlist_good(:).pred_nov_vs_fam]<0 & ~isnan([Neuronlist_good(:).rewardvalueindex_precue]));
else
    error('invalide Include type');
end



%% a mega figure summarizing the info task's finding.
figure;
% show that 3 reward related indices's measure is above chance
nsubplot(169,169, plotplacesety{1}, plotplacesetx{1}); set(gca,'ticklength',4*get(gca,'ticklength'))
Pvaluenames = {'rewardvalueindexP_precue', 'RewInfoAnticipIndexP_split'};
Pvaluelabel = {'Reward value', 'Infoanti'};

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


% a bar plot summarizing the correlation of correlation of Novelty index with the indices in NFL task and
% info tasks.
corr_permutexy = {[4,1], [4,2], [4,5],[4,6], [1,5], [1,6], [2,5], [2,6]};
corr_label = {'Sensory surprise', 'Recency', 'Reward value', 'Infoanti'...
    'Reward value', 'Infoanti', 'Reward value', 'Infoanti'};

corrvalue = zeros(size(corr_permutexy));
corrp = zeros(size(corr_permutexy));
corrstd = zeros(size(corr_permutexy));
corr_n = zeros(size(corr_permutexy));
for xyw = 1: numel(corr_permutexy)
    
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
if strcmpi(Include_type, 'Noveltyexcited')
    y_limit = [-0.1,0.15];
elseif strcmpi(Include_type, 'Noveltyinhibited')
    %y_limit = [-0.5,0.5];
    y_limit = [-0.1,0.15];
end
%plots
nsubplot(169,169, plotplacesety{2}, plotplacesetx{2}); set(gca,'ticklength',4*get(gca,'ticklength'))
bar([1:4],corrvalue(1:4), 'k');
errorbar([1:4],corrvalue(1:4), corrstd(1:4), '.');
ylabel('correlation with novelty index');
ylim(y_limit);
set(gca, 'xtick', [1:4], 'xticklabel', corr_label(1:4), 'XTickLabelRotation', -45);

nsubplot(169,169, plotplacesety{3}, plotplacesetx{3}); set(gca,'ticklength',4*get(gca,'ticklength'))
bar([1:2],corrvalue(4+(1:2)), 'k');
errorbar([1:2],corrvalue(4+(1:2)), corrstd(4+(1:2)), '.');
ylabel('correlation with sensory surprise index');
ylim(y_limit);
set(gca, 'xtick', [1:2], 'xticklabel', corr_label(4+(1:2)), 'XTickLabelRotation', -45);
pbaspect([3,5,1])

nsubplot(169,169, plotplacesety{4}, plotplacesetx{4}); set(gca,'ticklength',4*get(gca,'ticklength'))
bar([1:2],corrvalue(6+(1:2)), 'k');
errorbar([1:2],corrvalue(6+(1:2)), corrstd(6+(1:2)), '.');
ylabel('correlation with recency index');
ylim(y_limit);
set(gca, 'xtick', [1:2], 'xticklabel', corr_label(6+(1:2)), 'XTickLabelRotation', -45);
text(2,-0.08, ['n = ' mat2str(corr_n(1))]);
pbaspect([3,5,1])

% test to compare whether two rhos are significantly different.
nsubplot(169,169, plotplacesety{5}, plotplacesetx{5}); set(gca,'ticklength',4*get(gca,'ticklength'))

Corrset1 = {[4,1],[4,2]};
Corrset2 = {[4,5],[4,6]};
xtext = 0;
ytext = 0;
for ii=1:numel(Corrset1)
    for jj = 1:numel(Corrset2)
        
        xaxis_ind1 = [Neuronlist_good(Include_neurons).(indices_for_plot{Corrset1{ii}(1)})];
        yaxis_ind1 = [Neuronlist_good(Include_neurons).(indices_for_plot{Corrset1{ii}(2)})];
        xaxis_ind2 = [Neuronlist_good(Include_neurons).(indices_for_plot{Corrset2{jj}(1)})];
        yaxis_ind2 = [Neuronlist_good(Include_neurons).(indices_for_plot{Corrset2{jj}(2)})];
        
        
        % exclude the nan index
        excludeind = isnan(xaxis_ind1) | isnan(yaxis_ind1) | isnan(xaxis_ind2) | isnan(yaxis_ind2);
        xaxis_ind1(excludeind) = [];
        yaxis_ind1(excludeind) = [];
        xaxis_ind2(excludeind) = [];
        yaxis_ind2(excludeind) = [];
        p = Test_two_corr_difference_fun(xaxis_ind1, yaxis_ind1, xaxis_ind2, yaxis_ind2, shuffling_num);
        
        text(xtext, ytext, sprintf([axislabel_for_plot{Corrset1{ii}(1)} '+' axislabel_for_plot{Corrset1{ii}(2)} ' vs ' axislabel_for_plot{Corrset2{jj}(2)} 'p = %.5f']...
            , p), 'fontsize', 7);
        ytext = ytext+1;
    end
end
ylim([-1 ytext]);

nsubplot(169,169, plotplacesety{6}, plotplacesetx{6});
xlim([1,10]);
ylim([0,numel(corr_permutexy)+1]);
for xyw = 1: numel(corr_permutexy)
    text(1, xyw, [axislabel_for_plot{corr_permutexy{xyw}(1)} 'vs.' axislabel_for_plot{corr_permutexy{xyw}(2)} ', rho=' mat2str(corrvalue(xyw),4)...
        ', p = ' mat2str(corrp(xyw),4)]);
end



set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation



