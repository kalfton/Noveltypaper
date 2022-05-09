function Region_clustering_func(Neuronlist_all, neuron_counts, segmentmap, pars)
experimentConditions = pars.experimentConditions(1:3);
plotpath = pars.plotpath;
segmentNames = segmentmap(:,1);
allRegionsCells = neuron_counts.allRegionsCells;
significantCounts = neuron_counts.significantCounts;
significantPerc = neuron_counts.significantPerc;
pouts = neuron_counts.pouts;

%%% Matrices which stores the p value of pairwise comparison of regions
figure;
axplacey = {1:10, 11:20, 11:20, 1:10};
axplacex = {1:8, 11:18, 1:8, 11: 18};
saved_matrices = cell(3,1);
for condInd = 1:numel(experimentConditions)

ax = nsubplot(20,20,axplacey{condInd},axplacex{condInd});

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
color_star = [1 0 0];
for ii = 1:numel(segmentNames)
    for jj = 1:numel(segmentNames)
        %text(ii, jj, num2str(comparisonmatrix_p(ii, jj),2));
        if comparisonmatrix_p(ii, jj) <= pars.sigTres_low
            text(ii, jj, '***','Color',color_star, 'Fontsize', 5);
        elseif comparisonmatrix_p(ii, jj) <= pars.sigTres_mid
            text(ii, jj, '**','Color',color_star, 'Fontsize', 5);
        elseif comparisonmatrix_p(ii, jj) <= pars.sigTres_high
            text(ii, jj, '*','Color',color_star, 'Fontsize', 5);
        end
    end
end
saved_matrices{condInd} = comparisonmatrix;

title(pars.legend{condInd});
colorbar;

end

% statistics on the 3 matrices
nsubplot(20,20,axplacey{4},axplacex{4});
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
print(gcf,'-dpdf', '-painters',fullfile(plotpath, ['brain_area_pairwise_comparison_matrices_new' '.pdf']));

%% Cluster the regions by the percentage of neurons.
conditionnames = fieldnames(significantPerc.selective);
X_data = [];

for condind = 1:numel(conditionnames)
    X_data = [X_data, [significantPerc.selective.(conditionnames{condind})]'];
end

%% save the novelty/sensory surprise/recency indice of each neuron by brain areas
cells_areas = [Neuronlist_all.(pars.anatomyCase)];
areas_indices = cell(numel(segmentNames),1);
areas_indices_mean = zeros(numel(segmentNames),3);
for segmentInd = 1:numel(segmentNames)
    currentCells = ismember(cells_areas, segmentmap{segmentInd,2});
    for jj = 1:3 %novelty/sensory surprise/recency
        areas_indices{segmentInd}(:,end+1) = abs([Neuronlist_all(currentCells).(pars.experimentindices{jj})]);
        areas_indices_mean(segmentInd, jj) = mean(abs([Neuronlist_all(currentCells).(pars.experimentindices{jj})]));
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


set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);
print(gcf,'-dpdf', '-painters',fullfile(plotpath, ['brain_area_clustering' '.pdf']));

end