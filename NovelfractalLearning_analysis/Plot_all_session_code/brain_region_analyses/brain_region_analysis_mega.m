clear pars
pars.StatisticalThreshold = 0.01;
pars.sigThres_learning = 0.05;
pars.sigTres_low = 0.001;
pars.sigTres_mid = 0.01;
pars.sigTres_high = 0.05;
pars.experimentConditions = {'P_pred_nov_vs_fam','P_pred_vs_unpred_fam_perm', 'P_recency_ind_match_pos', 'P_violation_ind_perm'};
pars.Indexname = {'Novelty_ind', 'Surprise_ind', 'Recency_ind', 'Violation_ind'};
pars.counting_pairs = {[1], [2], [3], [4]};
pars.experimentindices = {'pred_nov_vs_fam', 'pred_vs_unpred_fam', 'recency_ind_match_pos', 'violation_ind'};
pars.experimentConditions_color = {'r','b','g', [0.8,0.2,0.8]};
pars.legend = {'Novelty', 'Sensory surprise', 'Recency', 'Violation'};
pars.cellThreshold = 90; % an area has to have minimum number of neurons to be included in the analysis 
pars.exclude_region = {'Claustrum', 'Zona Incerta'};
pars.anatomyCase = 'regionIndex';
pars.celltype = {'positive', 'negative','selective'};
pars.second_color = {'b','g',[0.8,0.2,0.8]};
pars.plotpath = plotpath;

if exist('shuffling_num', 'var')
    pars.shuffling_num = shuffling_num;
else
    pars.shuffling_num = 10000; % default shuffling number
end

experimentConditions = pars.experimentConditions;


%% count the percentage of neurons and raw numbers in each brain area.
[neuron_counts,segmentmap] = neuron_counts_region_V2(segmentmap, Neuronlist_all, pars);
%%
% allRegionsCells = neuron_counts.allRegionsCells;
% significantCounts = neuron_counts.significantCounts;
% significantPerc = neuron_counts.significantPerc;
% pouts = neuron_counts.pouts;
%plots

% sort the area by percentage of selective neuron
[~,order] = sort([neuron_counts.significantPerc.selective(1:end).P_pred_nov_vs_fam]);
%order = [numel(order),order]; % the first one should be alway not assigned.
neuron_counts.allRegionsCells = neuron_counts.allRegionsCells(order,1);
segmentmap = segmentmap(order,:);
celltype = pars.celltype;
for celltypeInd = 1:numel(celltype)
    neuron_counts.significantCounts.(celltype{celltypeInd}) = neuron_counts.significantCounts.(celltype{celltypeInd})(order);
    neuron_counts.significantPerc.(celltype{celltypeInd}) = neuron_counts.significantPerc.(celltype{celltypeInd})(order);
    neuron_counts.pouts.(celltype{celltypeInd}) = neuron_counts.pouts.(celltype{celltypeInd})(order);
end

% plot the percentage of neuron in each brain area and its correlation
region_perc_and_corr_func(Neuronlist_all, neuron_counts, segmentmap, pars);

% plot the comparison matrix and clustering analysis
Region_clustering_func(Neuronlist_all, neuron_counts, segmentmap, pars);

% plot the learning and forgetting by brain region analysis
learning_forgetting_region_func(Neuronlist_all, neuron_counts, segmentmap, pars);





