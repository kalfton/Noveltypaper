%% This script is for inter-region comparison
directories.resultsDir = plotpath;
monkeyNames = ["Lemmy","Slayer","Combined"];

pars.sigThres = 0.01;
pars.tableColor_low = 0.001;
pars.tableColor_mid = 0.01;
pars.tableColor_high = 0.05;
pars.experimentConditions = {'P_pred_nov_vs_fam', 'P_pred_vs_unpred_fam_perm', 'P_recency_ind_match_pos'};
pars.counting_pairs = {[1], [2], [3], [1,2], [2,3], [3,1], [1,2,3]};
correlation_pairs = {[1,2], [1,3]};
pars.experimentindices = {'pred_nov_vs_fam', 'pred_vs_unpred_fam', 'recency_ind_match_pos', 'violation_ind'};
pars.experimentConditions_color = {'r','b','g', [0.8,0.2,0.8]};
pars.legend = {'Novelty', 'Sensory surprise', 'Recency', 'Violation'};
pars.cellThreshold = 90;
pars.exclude_region = {'Claustrum', 'Zona Incerta'};
pars.include_region = {'AVMTE', '9/46V', 'Basal forebrain', 'Amygdala',...
    'Posterior Medial Temporal Cortex', '45B', 'OFC', 'Striatum'...
    'Anterior Entorhinal cortex', 'Globus Pallidus', 'Hippocampus'...
    'DIP', '8AD', 'Posterior Entorhinal cortex', 'Thalamus', '8A', '8B',...
    '6DC', '6DR', 'ACC', '4'};
pars.celltype = {'positive', 'negative','selective'};

monkeyName = 'Combined'; % 'Slayer','Lemmy'

pars.StatisticalThreshold=0.01;
pars.celltype = {'positive', 'negative','selective'};
pars.applyCellThreshold = 1;
pars.anatomyCase = 'regionIndex';
pars.normalize = true;

createneuronlist = true;
barinterv = 5;
%shuffling_num = 10000;

%parameters for correlation plot
Noveltytypes = {'Nov_excited', 'Nov_inhibited', 'Novel_selective','All_included'};
Indexname = {'Novelty_ind', 'Surprise_ind', 'Recency_ind', 'Violation_ind'};
correlation_legend = {'Novelty & Surprise', 'Novelty & Recency', 'Novelty & violation', 'Surprise & Recency'};
pars.correlation_color = {'b','g',[0.8,0.2,0.8]};
minmum_num = 10; % if an area has neuron number less than this, exclude it.


%%
try
    datastruct = Neuronlist_all;
catch
    datastruct = Neuronlist_good;
end

switch monkeyName
    case 'Lemmy'
        tempIndices = find([datastruct.MonkeyID]== 2);
        datastruct = datastruct(tempIndices);
    case 'Slayer'
        tempIndices = find([datastruct.MonkeyID]== 1);
        datastruct = datastruct(tempIndices);
    case 'Combined'
        
end

%clear Neuronlist_all;
experimentConditions = pars.experimentConditions;
StatisticalThreshold = pars.StatisticalThreshold;


%% count the percentage of neurons and raw numbers in each brain area.
[neuron_counts,segmentmap] = neuron_counts_region_V2(segmentmap, datastruct, pars);
segmentNames = segmentmap(:,1);
%%
allRegionsCells = neuron_counts.allRegionsCells;
significantCounts = neuron_counts.significantCounts;
significantPerc = neuron_counts.significantPerc;
pouts = neuron_counts.pouts;
%plots

% sort the area by percentage of selective neuron
[~,order] = sort([significantPerc.selective(1:end).P_pred_nov_vs_fam]);
%order = [numel(order),order]; % the first one should be alway not assigned.
allRegionsCells = allRegionsCells(order,1);
segmentNames = segmentNames(order);
segmentmap = segmentmap(order,:);
celltype = pars.celltype;
for celltypeInd = 1:numel(celltype)
    significantCounts.(celltype{celltypeInd}) = significantCounts.(celltype{celltypeInd})(order);
    significantPerc.(celltype{celltypeInd}) = significantPerc.(celltype{celltypeInd})(order);
    pouts.(celltype{celltypeInd}) = pouts.(celltype{celltypeInd})(order);
end


%% plot start here

%%% Matrices which stores the p value of pairwise comparison of regions
figure;
axplacey = {10:100, 110:200, 110:200, 10:100};
axplacex = {10:80, 110:180, 10:80, 110: 200};
saved_matrices = cell(3,1);
for condInd = 1:numel(pars.experimentConditions)

ax = nsubplot(210,210,axplacey{condInd},axplacex{condInd});

comparisonmatrix = zeros(numel(segmentNames));
comparisonmatrix_p = zeros(numel(segmentNames));
for ii = 1:numel(segmentNames)
    for jj = 1:numel(segmentNames)
        table_x = [significantCounts.selective(ii).(experimentConditions{condInd}), significantCounts.selective(jj).(experimentConditions{condInd})
            allRegionsCells(ii)-significantCounts.selective(ii).(experimentConditions{condInd}), allRegionsCells(jj)-significantCounts.selective(jj).(experimentConditions{condInd})];
        [~, comparisonmatrix_p(ii, jj)] = fishertest(table_x);
        comparisonmatrix(ii,jj) = abs(significantPerc.selective(ii).(experimentConditions{condInd}) - significantPerc.selective(jj).(experimentConditions{condInd}));
    end
end
colormap(bone);
imagesc(comparisonmatrix);

ax.XTick = 1:numel(segmentNames);
ax.YTick = 1:numel(segmentNames);
ax.XTickLabel = segmentNames;
ax.YTickLabel = segmentNames;
ax.XTickLabelRotation = -45;

% label the matrix with p value

for ii = 1:numel(segmentNames)
    for jj = 1:numel(segmentNames)
        %text(ii, jj, num2str(comparisonmatrix_p(ii, jj),2));
        if comparisonmatrix_p(ii, jj) <= pars.tableColor_low
            text(ii, jj, '***','Color',color_star, 'Fontsize', 5);
        elseif comparisonmatrix_p(ii, jj) <= pars.tableColor_mid
            text(ii, jj, '**','Color',color_star, 'Fontsize', 5);
        elseif comparisonmatrix_p(ii, jj) <= pars.tableColor_high
            text(ii, jj, '*','Color',color_star, 'Fontsize', 5);
        end
    end
end
saved_matrices{condInd} = comparisonmatrix;

title(pars.legend{condInd});
colorbar;

end

% statistics on the 3 matrices
nsubplot(210,210,axplacey{4},axplacex{4});
matrixpairs = {[1,2], [2,3], [1,3]};
for ii = 1:numel(matrixpairs)
    data1 = saved_matrices{matrixpairs{ii}(1)}(:);
    data2 = saved_matrices{matrixpairs{ii}(2)}(:);
    [rho, p] = corr(data1, data2, 'Type', 'Spearman');
    text(0,ii, [pars.legend{matrixpairs{ii}(1)} ' vs ' pars.legend{matrixpairs{ii}(2)} ', rho =' mat2str(rho,4) ', p=' mat2str(p,4)]);
end
ylim([-1, 10])
axis off;


set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
print(gcf,'-dpdf', '-painters',fullfile(directories.resultsDir, ['brain_area_pairwise_comparison_matrices_new' '.pdf']));

%% Cluster the regions by the percentage of neurons.
conditionnames = fieldnames(significantPerc.selective);
X_data = [];

for condind = 1:numel(conditionnames)
    X_data = [X_data, [significantPerc.selective.(conditionnames{condind})]'];
end

%% save the novelty/sensory surprise/recency indice of each neuron by brain areas
cells_areas = [datastruct.(pars.anatomyCase)];
areas_indices = cell(numel(segmentNames),1);
areas_indices_mean = zeros(numel(segmentNames),3);
for segmentInd = 1:numel(segmentNames)
    currentCells = ismember(cells_areas, segmentmap{segmentInd,2});
    for jj = 1:3 %novelty/sensory surprise/recency
        areas_indices{segmentInd}(:,end+1) = abs([datastruct(currentCells).(pars.experimentindices{jj})]);
        areas_indices_mean(segmentInd, jj) = mean(abs([datastruct(currentCells).(pars.experimentindices{jj})]));
    end
    clear currentCells;
end

distances = zeros(1, numel(segmentNames)*(numel(segmentNames)-1)/2);
zz = 1;
for ii = 1:numel(segmentNames)-1
    for jj = ii+1:numel(segmentNames)
        dist_matrix = pdist2(areas_indices{ii}, areas_indices{jj}, 'euclidean');
        distances(zz) = mean(dist_matrix(:));
        zz = zz+1;
    end
end
% normalize areas_indices_mean
areas_indices_mean_n = areas_indices_mean./areas_indices_mean(:,1);

distances_mean = pdist(areas_indices_mean_n);
distances_percent = pdist(X_data);

hierachy_cluster = linkage(distances_mean, 'average');
figure;
graph1 = subplot(4,2,[1,3]);
[plotobj, ~, dendr_perm] = dendrogram(hierachy_cluster, 'Labels', segmentNames);
set(gca,'XTickLabelRotation',-25);
xl = xlim;

plotInd = 0;
barInd = zeros(1,numel(experimentConditions));
clear plot_id;

% normalize the mean of the indices
if pars.normalize
    areas_indices_mean_n2 = (areas_indices_mean-min(areas_indices_mean,[],1))./(max(areas_indices_mean,[],1)-min(areas_indices_mean,[],1));
else
    areas_indices_mean_n2 = areas_indices_mean;
end
celltype = pars.celltype;
graph_c = subplot(4,2,5);
for condInd = 1:numel(experimentConditions)
    plotInd = plotInd + 1;
    w = 0.05;
    color = [0 0.7 0.7];
    xdata = areas_indices_mean_n2(dendr_perm,condInd);
    percCoord = xdata;
    barInd(condInd) = plotInd;
    color_star = [1 0 0];
    barposition = (barinterv)*(1:numel(xdata))+numel(experimentConditions)-condInd;
    %plot_id(plotInd) = barh(barposition, perc,w,'FaceColor',pars.experimentConditions_color{condInd});
    %plot the trending line
    %barposition = 1:numel(xdata);
    plot_id(plotInd) = plot(barposition, xdata, 'Linewidth', 1, 'Color',pars.experimentConditions_color{condInd});
    hold on;
end

xlim(xl*barinterv);
graph_c.Position(1) = graph1.Position(1);
graph_c.Position(3) = graph1.Position(3);
ylim([0,1]);

ylabel('Normalized averaged indices');

set(gca,'Xtick', barposition, 'XtickLabel',segmentNames(dendr_perm),'FontAngle','italic');
set(gca,'XTickLabelRotation',-25);
set(gca,'Ytick', [0:0.2:1]);

% another plot which plot the none normalized averaged roc value.
graph_c = subplot(4,2,7);
celltype = pars.celltype;
for condInd = 1:numel(experimentConditions)
    plotInd = plotInd + 1;
    w = 0.05;
    color = [0 0.7 0.7];
    xdata = areas_indices_mean(dendr_perm,condInd);
    percCoord = xdata;
    barInd(condInd) = plotInd;
    color_star = [1 0 0];
    barposition = (barinterv)*(1:numel(xdata))+numel(experimentConditions)-condInd;
    %plot the trending line
    plot_id(plotInd) = plot(barposition, xdata, 'Linewidth', 1, 'Color',pars.experimentConditions_color{condInd});
    hold on;
end

xlim(xl*barinterv);
graph_c.Position(1) = graph1.Position(1);
graph_c.Position(3) = graph1.Position(3);

ylim([0,0.4]); 

legend(plot_id(barInd),pars.legend,'Location','northwest');
ylabel('Averaged indices');
xlabel('Segment Names');

set(gca,'Xtick', barposition, 'XtickLabel',segmentNames(dendr_perm),'FontAngle','italic');
set(gca,'XTickLabelRotation',-25);
set(gca,'Ytick', [0:0.2:1]);


set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);
print(gcf,'-dpdf', '-painters',fullfile(directories.resultsDir, ['brain_area_normalized_ROC' '.pdf']));







