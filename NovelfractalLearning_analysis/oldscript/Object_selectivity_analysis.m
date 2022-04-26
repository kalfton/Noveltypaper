% analysis of the correlation of timescale and the object selectivity index

% object selectivity index comparison
figure;

indice_1_name = {'obj_select_ind_novel_start', 'obj_select_ind_recent', 'obj_select_ind_surprise'};
indice_2_name = {'obj_select_ind_novel_end', 'obj_select_ind_nonrecent', 'obj_select_ind_nonsurprise'};
indice_1_alias = {'novel start', 'recent', 'surprise'};
indice_2_alias = {'novel end', 'nonrecent', 'nonsurprise'};
ploty = {10:50, 60:100, 110:150};
plotx = {10:70, 10:70, 10:70};
titlename = {'Novelty', 'Recency', 'Sensory surprise'};

include_neurons = true(size(Neuronlist_good));%[Neuronlist_good.P_pred_nov_vs_fam]'<StatisticalThreshold;% & [Neuronlist_good(:).pred_nov_vs_fam]'>0

for ii = 1:numel(titlename)
obj_ind1 = [Neuronlist_good(include_neurons).(indice_1_name{ii})];
obj_ind2 = [Neuronlist_good(include_neurons).(indice_2_name{ii})];

nanind = isnan(obj_ind1) | isnan(obj_ind2);
obj_ind1(nanind)=[];
obj_ind2(nanind)=[];
% sanitiy check of the distribution of the object selectivity indices
nsubplot(169, 169, ploty{ii}, plotx{ii}+80);
histogram(obj_ind1,'FaceColor', 'r');hold on;
histogram(obj_ind2, 'FaceColor', 'b');
xlabel('object selectivity index');

nsubplot(169, 169, ploty{ii}, plotx{ii});
n = numel(obj_ind1);
ind_mean = [mean(obj_ind1), mean(obj_ind2)];
ind_std = [std(obj_ind1)/sqrt(n), std(obj_ind2)/sqrt(n)];

bar(1, ind_mean(1), 'r');
bar(2, ind_mean(2), 'b'); 
errorbar([1,2], ind_mean, ind_std, '.');
set(gca, 'xtick', [1,2], 'xticklabel', {indice_1_alias{ii},indice_2_alias{ii}});
ylabel('object selectivity index');
% satistical test
 p = signrank(obj_ind1, obj_ind2);
title([titlename{ii}, 'n = ' mat2str(numel(obj_ind1)) ', p = ' mat2str(p, 4)]);

end

set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
print(gcf,'-dpdf', '-painters',fullfile(plotpath,'object_selectivity analysis.pdf'));