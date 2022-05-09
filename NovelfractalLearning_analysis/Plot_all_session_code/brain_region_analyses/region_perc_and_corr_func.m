function region_perc_and_corr_func(Neuronlist_all, neuron_counts, segmentmap, pars)
StatisticalThreshold=pars.StatisticalThreshold;
shuffling_num = pars.shuffling_num;
Indexname = pars.Indexname;
celltype = pars.celltype;
experimentConditions = pars.experimentConditions;
plotpath = pars.plotpath;
correlation_pairs = {[1,2], [1,3], [1,4]};
correlation_legend = {'Novelty & Surprise', 'Novelty & Recency', 'Novelty & violation'};
minmum_num = 20; % if an area has neuron number less than this, exclude it.
barinterv = 5;
Noveltytypes = 'All_included'; %'Nov_excited', 'Nov_inhibited'

segmentNames = segmentmap(:,1);
allRegionsCells = neuron_counts.allRegionsCells;
significantCounts = neuron_counts.significantCounts;
significantPerc = neuron_counts.significantPerc;
pouts = neuron_counts.pouts;

%% Ploting start here.

figure;
nsubplot(21,21,6:20,1:10);
plotInd = 0;
barInd = zeros(1,numel(experimentConditions));
clear plot_id;

for celltypeInd = numel(celltype) %1:2%numel(celltype)
currentCellCase = celltype{celltypeInd};
for condInd = 1:numel(experimentConditions)
    plotInd = plotInd + 1;
    perc = 100 * [significantPerc.(currentCellCase)(:).(experimentConditions{condInd})];
    if strcmpi(currentCellCase, 'negative')
        perc = -perc;
    end
    percCoord = perc;
    barInd(condInd) = plotInd;
    barposition = (barinterv)*(1:numel(perc))+numel(experimentConditions)-condInd;
    %plot the trending line
    plot_id(plotInd) = plot(perc,barposition,  'Linewidth', 1, 'Color',pars.experimentConditions_color{condInd});
    hold on;
    if strcmpi(currentCellCase, 'selective')
        xlim([0 60]);
    else
        xlim([-60 60]);
    end
    xl = xlim ; xl = xl(end);
    adj = xl/30;
    if strcmpi(currentCellCase, 'negative')
        adj = -adj;
    end
    
    if condInd == 1 && contains( currentCellCase, {'positive', 'selective'})
        for segmentInd = 1:numel(segmentNames)
            text(percCoord(segmentInd) + adj, barposition(segmentInd)-floor(numel(experimentConditions)/2), num2str(significantCounts.(celltype{celltypeInd})(segmentInd).(experimentConditions{condInd})));
            %  text(percCoord(segmentInd) + 2*adj, segmentInd,num2str(significantCounts(segmentInd).('allSig')));
            
            text(percCoord(segmentInd) + 3*adj, barposition(segmentInd)-floor(numel(experimentConditions)/2), num2str(allRegionsCells(segmentInd)));
        end
    end   

end

% an additional analysis to ask whether neurons of novelty and other
% indices are correlated or not
perc_array_novelty = [significantPerc.(currentCellCase)(:).(experimentConditions{1})];
perc_array_surprise = [significantPerc.(currentCellCase)(:).(experimentConditions{2})];
perc_array_recency = [significantPerc.(currentCellCase)(:).(experimentConditions{3})];

if strcmpi(currentCellCase, 'negative')
    textx = -40;
elseif strcmpi(currentCellCase, 'positive')
    textx = 40;
else
    textx = 40;
end

[rho, p] = corr(perc_array_novelty', perc_array_surprise', 'Type', 'Spearman');
text(textx, barposition(end), sprintf('novel& surprise, rho = %.4f, p = %.4f', rho,p));
[rho, p] = corr(perc_array_novelty', perc_array_recency', 'Type', 'Spearman');
text(textx, barposition(end-1), sprintf('novel& recency, rho = %.4f, p = %.4f', rho,p));

% another statistical test
n_divid = ceil(numel(perc_array_novelty)/3);
p = ranksum(perc_array_surprise(1:n_divid), perc_array_surprise(end-n_divid+1:end));
text(textx, barposition(end-2), sprintf('surprise ranksum, p = %.4f',p));
p = ranksum(perc_array_recency(1:n_divid), perc_array_recency(end-n_divid+1:end));
text(textx, barposition(end-3), sprintf('recency ranksum, p = %.4f',p));

%yet another statistical test
if strcmpi(currentCellCase, 'selective')
    p_thrh = 0.01;
else
    p_thrh = 0.01/2;
end
% total number of brain areas
n_area = numel(significantCounts.(celltype{celltypeInd}));

% calculate the top 4 brain areas and the lowest 4 brain area using binomial test
% lowest 4:
s = 0;
n = 0;
for ii = 1:4
s = s + significantCounts.(celltype{celltypeInd})(ii).('P_recency_ind_match_pos');
n = n + allRegionsCells(ii);
end
s_lowest = s;
n_lowest = n;
p_recency_lowest = myBinomTest(s,n,p_thrh,'Greater');
text(textx, barposition(end-5), sprintf('recency, lowest 4 area, p = %.4f',p_recency_lowest));

% top 4:
s = 0;
n = 0;
for ii = n_area-4+1:n_area
s = s + significantCounts.(celltype{celltypeInd})(ii).('P_recency_ind_match_pos');
n = n + allRegionsCells(ii);
end
s_top = s;
n_top = n;
p_recency_top = myBinomTest(s,n,p_thrh,'Greater');
text(textx, barposition(end-6), sprintf('recency, highest 4 area, p = %.4f',p_recency_top));

p_ratio_compare = permutation_exact_test_of_proportions(s_top, n_top, s_lowest, n_lowest);
text(textx, barposition(end-7), sprintf('recency, highest vs lowest, p = %.4f',p_ratio_compare));

% lowest 4:
s = 0;
n = 0;
for ii = 1:4
s = s + significantCounts.(celltype{celltypeInd})(ii).('P_pred_vs_unpred_fam_perm');
n = n + allRegionsCells(ii);
end
s_lowest = s;
n_lowest = n;
p_surprise_lowest = myBinomTest(s,n,p_thrh,'Greater');
text(textx, barposition(end-8), sprintf('surprise, lowest 4 area, p = %.4f',p_surprise_lowest));

% top 4:
s = 0;
n = 0;
for ii = n_area-4+1:n_area
s = s + significantCounts.(celltype{celltypeInd})(ii).('P_pred_vs_unpred_fam_perm');
n = n + allRegionsCells(ii);
end
s_top = s;
n_top = n;
p_surprise_top = myBinomTest(s,n,p_thrh,'Greater');
text(textx, barposition(end-9), sprintf('surprise, highest 4 area, p = %.4f',p_surprise_top));

p_ratio_compare = permutation_exact_test_of_proportions(s_top, n_top, s_lowest, n_lowest);
text(textx, barposition(end-10), sprintf('surprise, highest vs lowest, p = %.4f',p_ratio_compare));

end

legend(plot_id(barInd),pars.legend,'Position',[0.4 0.9 0.05 0.025]);
xlabel('% of cells');
ylabel('Segment Names');

str2title = {[' Single Units Original Areas '] ; ['# of cells']};
str2title{end + 1} = ['Including only regions with cell # >',num2str(pars.cellThreshold)];
title(str2title);

yticks(1:numel(segmentNames));
set(gca,'Ytick', (barinterv)*(1:numel(segmentNames))+floor(numel(experimentConditions)/2), 'YtickLabel',segmentNames,'FontAngle','italic');
set(gca,'Xtick', [-100:20:100], 'xtickLabel',[100:-20:0,20:20:100]);



%%
clear included_Ind included_Ind_pos included_Ind_neg
if contains(Noveltytypes, 'Nov_excited')
    included_Ind = find([Neuronlist_all.P_pred_nov_vs_fam] <= StatisticalThreshold & [Neuronlist_all.pred_nov_vs_fam] >0 );
elseif contains(Noveltytypes, 'Nov_inhibited')
    included_Ind = find([Neuronlist_all.P_pred_nov_vs_fam] <= StatisticalThreshold & [Neuronlist_all.pred_nov_vs_fam] <0 );
elseif contains(Noveltytypes, 'All_included')
    included_Ind = 1:numel(Neuronlist_all);%find([Neuronlist_all.P_pred_nov_vs_fam] <= inf );
else
    error('Wrong Novelty type');
end

clear Indices_struct
for ii = 1: numel(Indexname)
    Indices_struct.(Indexname{ii}) = [Neuronlist_all(included_Ind).(pars.experimentindices{ii})];
end


%% count the # of cells with current type in each area
clear allRegionsCells
cells_areas = [Neuronlist_all(included_Ind).(pars.anatomyCase)];

for segmentInd = 1:numel(segmentNames)
    currentSegmentNo = segmentmap{segmentInd,2};
    currentCells = ismember(cells_areas, currentSegmentNo);
    allRegionsCells(segmentInd,1) = numel(find(currentCells == 1));
    clear currentCells;
end

%% caluculate correlations in each area
clear Correlations P_corr Correlation_bounds
cells_areas = [Neuronlist_all(included_Ind).(pars.anatomyCase)];
for condInd = 1:numel(correlation_pairs)
    for segmentInd = 1:numel(segmentNames) % calculate correlation in each brain area
        currentSegmentNo = segmentmap{segmentInd,2};
        
        Index1 = Indices_struct.(Indexname{correlation_pairs{condInd}(1)});
        Index2 = Indices_struct.(Indexname{correlation_pairs{condInd}(2)});
        Index1 = Index1(ismember(cells_areas, currentSegmentNo));
        Index2 = Index2(ismember(cells_areas, currentSegmentNo));
        % delete the neurons that has nan in either of the indices
        nanindex = isnan(Index1) | isnan(Index2);
        Index1(nanindex) = [];
        Index2(nanindex) = [];
        
        if numel(Index1)<minmum_num 
            rho = nan;
            p_val = nan;
            rho_upperbound = nan;
            rho_lowerbound = nan;
            rho_shuffling = zeros(shuffling_num,1);
        else
             [rho, p_val] = corr(Index1',Index2', 'Type', 'Spearman');
            % bootstrapping to get confidence interval
            rho_shuffling = zeros(shuffling_num,1);
            
            for ii = 1:shuffling_num
                shuffle_Ind = randi(numel(Index1),(size(Index1)));
                Index1_shuffle = Index1(shuffle_Ind);
                Index2_shuffle = Index2(shuffle_Ind);
                rho_shuffling(ii) = corr(Index1_shuffle',Index2_shuffle', 'Type', 'Spearman');
            end
            rho_shuffling = sort(rho_shuffling);
            rho_upperbound = rho_shuffling(shuffling_num-floor(0.05*shuffling_num));
            rho_lowerbound = rho_shuffling(floor(0.05*shuffling_num));
            
        end
        
        Correlations.(Noveltytypes)(segmentInd,condInd) = rho;
        P_corr.(Noveltytypes)(segmentInd,condInd) = p_val;
        % use std of the shuffling data.
        Correlation_bounds.(Noveltytypes).upperbound(segmentInd,condInd) = rho+std(rho_shuffling);
        Correlation_bounds.(Noveltytypes).lowerbound(segmentInd,condInd) = rho-std(rho_shuffling);
    end
end

subplotposition_x = {11:13, 14:16,17:19};
subplotposition_y = {1:5, 6:20};

plotInd = 0;
barInd = zeros(1,numel(correlation_pairs));
clear plot_id

for condInd = 1:numel(correlation_pairs)
    nsubplot(21,21,subplotposition_y{2}, subplotposition_x{condInd});
    plotInd = plotInd + 1;
    perc = [Correlations.(Noveltytypes)(:,condInd)];
    upperbound = [Correlation_bounds.(Noveltytypes).upperbound(:,condInd)]-perc;
    lowerbound = -([Correlation_bounds.(Noveltytypes).lowerbound(:,condInd)]-perc);
    barInd(condInd) = plotInd;
    barposition = (barinterv)*(1:numel(perc))+numel(correlation_pairs)-condInd;
    hold on;
    %error bar
    errorbar(perc, barposition,(upperbound+lowerbound)/2, -(upperbound+lowerbound)/2, 'horizontal', 'LineStyle','none', 'Marker', '.', 'color', pars.second_color{condInd});
    %plot the trending line
    plot_id(plotInd) = plot(perc,barposition,  'Linewidth', 1, 'Color',pars.second_color{condInd});
    xlim([-0.4 0.4]);
    xl = xlim ; xl = xl(end);
    xlabel('Correlation (rho)');
if condInd==1
    set(gca,'Ytick', (barinterv)*(1:numel(segmentNames))+floor(numel(correlation_pairs)/2), 'YtickLabel',segmentNames,'FontAngle','italic');
else
    set(gca,'Ytick', (barinterv)*(1:numel(segmentNames))+floor(numel(correlation_pairs)/2), 'YtickLabel',[]);
end
end

%% Plot an histogram on top of the correlation bar plot
for condInd = 1:numel(correlation_pairs)
    nsubplot(21,21,subplotposition_y{1}, subplotposition_x{condInd});
    perc = [Correlations.(Noveltytypes)(:,condInd)];
    edges = [-0.4:0.05:0.4];
    histogram(perc,edges,'facecolor', pars.second_color{condInd});
    hold on;
    xlim([-0.4 0.4]);
    ylim([0,15]);
    %mean of the histogram
    mean(perc);
    YL = get(gca, 'ylim');
    YL = YL(end);
    scatter(nanmean(perc),YL,80,pars.second_color{condInd},'v','filled')
    line([nanmean(perc) nanmean(perc)], [0 YL] , 'Color', pars.second_color{condInd},'LineWidth',1) ;
    
    p = signrank(perc,0);
    text(-0.3 ,YL-2*condInd, sprintf('p = %.5f',p), 'Color', pars.second_color{condInd});
    pbaspect([7,4,1])
end

legend(plot_id(barInd),correlation_legend(1:numel(barInd)),'Position',[0.9 0.9 0.05 0.025]);

set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);
print(fullfile(plotpath,['Main_neuron_distribution_and_correlation_across_area_plot.pdf']),'-dpdf');

end