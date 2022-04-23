%% This is for plotting the learning and forgetting by brain area plot
monkeyNames = ["Lemmy","Slayer","Combined"];

clear pars;
pars.sigThres = 0.01;
pars.sigThres_learning = 0.05;
pars.tableColor_low = 0.001;
pars.tableColor_mid = 0.01;
pars.tableColor_high = 0.05;
pars.experimentConditions = {'P_pred_nov_vs_fam','P_pred_vs_unpred_fam_perm', 'P_recency_ind_match_pos'}; %'P_violation_ind_perm' 'P_pred_vs_unpred_fam_perm', 'P_recency_ind_match_pos'
pars.counting_pairs = {[1], [2], [3]};
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
pars.separate_pos_neg = true;
pars.ifCombined = 0;
pars.applyCellThreshold = 1;
pars.anatomyCase = 'regionIndex';


anatomyCase = 'regionIndex';
createneuronlist = false;
StatisticalThreshold=0.01;
barinterv = 5;

%parameters for correlation plot
Noveltytypes = {'Nov_excited', 'Nov_inhibited', 'Novel_selective','All_included'};
Indexname = {'Novelty_ind', 'Surprise_ind', 'Recency_ind', 'Violation_ind'};
correlation_pairs = {[1,2], [1,3]};
correlation_legend = {'Novelty & Surprise', 'Novelty & Recency', 'Novelty & violation', 'Surprise & Recency'};
pars.correlation_color = {'b','g',[0.8,0.2,0.8]};
minmum_num = 20; % if an area has neuron number less than this, exclude it.


%% Ploting start here.
experimentConditions = pars.experimentConditions;
%load(fullfile(directories.variableDir,'segmentmap_Kaining.mat'));% load segmentNames

%% count the percentage of neurons and raw numbers in each brain area.
[neuron_counts,segmentmap] = neuron_counts_region_V2(segmentmap, Neuronlist_good, pars);
segmentNames = segmentmap(:,1);
%%
allRegionsCells = neuron_counts.allRegionsCells;
significantCounts = neuron_counts.significantCounts;
significantPerc = neuron_counts.significantPerc;
pouts = neuron_counts.pouts;

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


%plots

figure;
nsubplot(1,2,1:1,1:1);
plotInd = 0;
barInd = zeros(1,numel(experimentConditions));
clear plot_id;

for celltypeInd = 1:2%1:2%numel(celltype)
currentCellCase = celltype{celltypeInd};
for condInd = 1:numel(experimentConditions)
    plotInd = plotInd + 1;
    w = 0.05;
    color = [0 0.7 0.7];
    perc = 100 * [significantPerc.(currentCellCase)(:).(experimentConditions{condInd})];
    if strcmpi(currentCellCase, 'negative')
        perc = -perc;
    end
    percCoord = perc;
    barInd(condInd) = plotInd;
    color_star = [1 0 0];
    barposition = (barinterv)*(1:numel(perc))+numel(experimentConditions)-condInd;
    plot_id(plotInd) = barh(barposition, perc,w,'FaceColor',pars.experimentConditions_color{condInd});
    %plot the trending line
    plot(perc,barposition,  'Linewidth', 1, 'Color',pars.experimentConditions_color{condInd});
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
    %%text
    %% percentage acc. to all cells
    for segmentInd = 1:numel(segmentNames)
        % write the actual #
        % sig of binomial test
        if pouts.(currentCellCase)(segmentInd).(experimentConditions{condInd}) <= pars.tableColor_low
            text(perc(segmentInd) + adj/2, barposition(segmentInd) + w/3,'***','Color',color_star);
        elseif pouts.(currentCellCase)(segmentInd).(experimentConditions{condInd}) <= pars.tableColor_mid
            text(perc(segmentInd) + adj/2, barposition(segmentInd) + w/3,'**','Color',color_star);
        elseif pouts.(currentCellCase)(segmentInd).(experimentConditions{condInd}) <= pars.tableColor_high
            text(perc(segmentInd) + adj/2, barposition(segmentInd) + w/3,'*','Color',color_star);
        end
        %
    end
    %
    xt = xticks;
    xEnd = xt(end);
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
    textx = -50;
elseif strcmpi(currentCellCase, 'positive')
    textx = 50;
else
    textx = 50;
end

[rho, p] = corr(perc_array_novelty', perc_array_surprise', 'Type', 'Spearman');
text(textx, barposition(end), sprintf('novel& surprise, rho = %.3f, p = %.3f', rho,p));
[rho, p] = corr(perc_array_novelty', perc_array_recency', 'Type', 'Spearman');
text(textx, barposition(end-1), sprintf('novel& recency, rho = %.3f, p = %.3f', rho,p));

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

if excludemua
    str2title = {[' Single Units Original Areas '] ; ['# of cells']};
else
    str2title = {[' Multi+Single Units Original Areas '] ; ['# of cells']};
end
if pars.applyCellThreshold
    str2title{end + 1} = ['Including only regions with cell # >',num2str(pars.cellThreshold)];
else
    str2title{end + 1} = ['No cell # threshold applied'];
end
title(str2title);

yticks(1:numel(segmentNames));
set(gca,'Ytick', (barinterv)*(1:numel(segmentNames))+floor(numel(experimentConditions)/2), 'YtickLabel',segmentNames,'FontAngle','italic');
set(gca,'Xtick', [-100:20:100], 'xtickLabel',[100:-20:0,20:20:100]);




%% learning index and forgetting index plot
currenttype = Noveltytypes{1};
included_Ind = find([Neuronlist_good.learningforgetinganalysis]);
Indexname = {'withindaylearning_newindex', 'acrossdayforget_newindex'};
Indexlegend = {'learning index', 'forgetting index'};

clear Indices_struct 
for ii = 1: numel(Indexname)
    Indices_struct.(Indexname{ii}) = [Neuronlist_good(included_Ind).(Indexname{ii})];
end

%% count the # of cells with current type in each area
clear allRegionsCells
cells_areas = [Neuronlist_good(included_Ind).(anatomyCase)];

for segmentInd = 1:numel(segmentNames)
    currentSegmentNo = segmentmap{segmentInd,2};
    currentCells = ismember(cells_areas, currentSegmentNo);
    allRegionsCells(segmentInd,1) = numel(find(currentCells == 1));
    clear currentCells;
end

%% average the indices within each area

clear Indexmean P_val Indexstd
cells_areas = [Neuronlist_good(included_Ind).(anatomyCase)];
for condInd = 1:numel(Indexname)
    for segmentInd = 1:numel(segmentNames) % calculate correlation in each brain area
        currentSegmentNo = segmentmap{segmentInd,2};
        
        current_index = Indices_struct.(Indexname{condInd});
        current_index = current_index(ismember(cells_areas, currentSegmentNo));
        %current_index = current_index(~isnan(current_index));
        
        if numel(current_index)<minmum_num 
            meanind = nan;
            p_val = nan;
            stdind = nan;
        else
             meanind = mean(current_index);
             stdind = std(current_index)/sqrt(numel(current_index));
             p_val = signrank(current_index);
            
        end
        
        Indexmean(segmentInd,condInd) = meanind;
        P_val(segmentInd,condInd) = p_val;
        Indexstd(segmentInd,condInd) = stdind;
        
    end
end

%% plot

nsubplot(3,2,1:2,2:2);
plotInd = 0;
barInd = zeros(1,numel(Indexname));
clear plot_id

% only include areas that have enough neurons
included_brainarea = all(~isnan(Indexmean), 2);
P_val = P_val(included_brainarea,:);
Indexmean = Indexmean(included_brainarea,:);
Indexstd = Indexstd(included_brainarea,:);
segmentNames = segmentNames(included_brainarea,:);
allRegionsCells = allRegionsCells(included_brainarea,:);

for condInd = 1:numel(Indexname)
    plotInd = plotInd + 1;
    w = 0.05;
    color = [0 0.7 0.7];
    
    perc = [Indexmean(:,condInd)];
    upperbound = Indexstd(:,condInd);
    lowerbound = Indexstd(:,condInd);
    percCoord = perc;
    barInd(condInd) = plotInd;
    color_star = [1 0 0];
    barposition = (barinterv)*(1:numel(perc))+numel(correlation_pairs)-condInd;
    plot_id(plotInd) = barh(barposition, perc,w,'FaceColor',pars.correlation_color{condInd});
    hold on;
    %error bar
    errorbar(perc, barposition,(upperbound+lowerbound)/2, -(upperbound+lowerbound)/2, 'horizontal', 'LineStyle','none', 'Marker', '.', 'color', pars.correlation_color{condInd});
    %plot the trending line
    plot(perc,barposition,  'Linewidth', 1, 'Color',pars.correlation_color{condInd});
    xlim([-0.05 0.1]);
    xl = xlim ; xl = xl(end);
    adj = xl/30;
    %% percentage acc. to all cells
    for segmentInd = 1:numel(segmentNames)
        % write the actual #
        % sig of binomial test
        if P_val(segmentInd,condInd) <= pars.tableColor_low
            text(perc(segmentInd) + adj/2, barposition(segmentInd) + w/3,'***','Color',color_star);
        elseif P_val(segmentInd,condInd) <= pars.tableColor_mid
            text(perc(segmentInd) + adj/2, barposition(segmentInd) + w/3,'**','Color',color_star);
        elseif P_val(segmentInd,condInd) <= pars.tableColor_high
            text(perc(segmentInd) + adj/2, barposition(segmentInd) + w/3,'*','Color',color_star);
        end
        %
        text(perc(segmentInd) + adj, barposition(segmentInd) + w/3 ,['p = ' mat2str(P_val(segmentInd,condInd),3)]);
    end
    %
    xt = xticks;
    xEnd = xt(end);
    if condInd == 1
        for segmentInd = 1:numel(segmentNames)
            text(percCoord(segmentInd) + 4*adj, barposition(segmentInd)-floor(numel(correlation_pairs)/2), num2str(allRegionsCells(segmentInd)));
        end
    end
end

textx = 0.05;
notnanlogical = ~isnan(Indexmean(:,1));
[rho, p] = corr(Indexmean(notnanlogical,1), Indexmean(notnanlogical,2), 'Type', 'Spearman');
text(textx, barposition(end), sprintf('learning & forgetting, rho = %.3f, p = %.3f', rho,p));

legend(plot_id(barInd),Indexlegend(1:numel(barInd)),'Position',[0.9 0.9 0.05 0.025]);
xlabel('Indices value');
ylabel('Segment Names');

if excludemua
    str2title = {[' Single Units Original Areas ']};
else
    str2title = {[' Multi+Single Units Original Areas '] };
end

str2title{end + 1} = ['Learning index and Forgetting index are averaged in areas with cell # >',num2str(minmum_num)];
str2title{end + 1} = ['Cell type: ',currenttype];


title(str2title, 'interpreter','none');



yticks(1:numel(segmentNames));
set(gca,'Ytick', (barinterv)*(1:numel(segmentNames))+floor(numel(correlation_pairs)/2), 'YtickLabel',segmentNames,'FontAngle','italic');
%set(gcf, 'Position', get(0, 'Screensize'));

set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);

nsubplot(3,2,3,2:2);
%set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);
scatter(Indexmean(:,1), Indexmean(:,2), 'filled')
xl = xlim ; xl = xl(end);
adj = xl/60;
for ii = 1:numel(segmentNames)
    text(Indexmean(ii,1)+adj, Indexmean(ii,2),segmentNames{ii}, 'Fontsize', 7);
end
xlabel('within day learning index');
ylabel('across day forgetting index');
xlim([-0.05, 0.1]);
ylim([-0.05, 0.1]);
plot(get(gca, 'xlim'),[0,0], 'k');
plot([0,0], get(gca, 'ylim'), 'k');
axis square;


% a bar pie version of the plot
%% count the # of cells with that is significant in learning, forgetting or both
clear allRegionsCells learning_only_Cells forgetting_only_Cells Both_Cells Neither_Cells
segmentNames = segmentmap(:,1);
cells_areas = [Neuronlist_good(included_Ind).(anatomyCase)];

for segmentInd = 1:numel(segmentNames)
    currentSegmentNo = segmentmap{segmentInd,2};
    currentCells = ismember(cells_areas, currentSegmentNo);
    allRegionsCells(segmentInd,1) = numel(find(currentCells == 1));
    learning_only_Cells(segmentInd,1) = nnz([Neuronlist_good(included_Ind).withindaylearning_newindex_p]<pars.sigThres_learning & [Neuronlist_good(included_Ind).acrossdayforget_newindex_p]>=pars.sigThres_learning & currentCells);
    forgetting_only_Cells(segmentInd,1) = nnz([Neuronlist_good(included_Ind).withindaylearning_newindex_p]>=pars.sigThres_learning & [Neuronlist_good(included_Ind).acrossdayforget_newindex_p]<pars.sigThres_learning & currentCells);
    Both_Cells(segmentInd,1) = nnz([Neuronlist_good(included_Ind).withindaylearning_newindex_p]<pars.sigThres_learning & [Neuronlist_good(included_Ind).acrossdayforget_newindex_p]<pars.sigThres_learning & currentCells);
    Neither_Cells(segmentInd,1) = nnz([Neuronlist_good(included_Ind).withindaylearning_newindex_p]>=pars.sigThres_learning & [Neuronlist_good(included_Ind).acrossdayforget_newindex_p]>=pars.sigThres_learning & currentCells);
    forgetting_pos_Cells(segmentInd,1) = nnz([Neuronlist_good(included_Ind).acrossdayforget_newindex]>=0 & [Neuronlist_good(included_Ind).acrossdayforget_newindex_p]<pars.sigThres_learning & currentCells);
    forgetting_neg_Cells(segmentInd,1) = nnz([Neuronlist_good(included_Ind).acrossdayforget_newindex]<0 & [Neuronlist_good(included_Ind).acrossdayforget_newindex_p]<pars.sigThres_learning & currentCells);
    learning_pos_Cells(segmentInd,1) = nnz([Neuronlist_good(included_Ind).withindaylearning_newindex]>=0 & [Neuronlist_good(included_Ind).withindaylearning_newindex_p]<pars.sigThres_learning & currentCells);
    learning_neg_Cells(segmentInd,1) = nnz([Neuronlist_good(included_Ind).withindaylearning_newindex]<0 & [Neuronlist_good(included_Ind).withindaylearning_newindex_p]<pars.sigThres_learning & currentCells);
    clear currentCells;
end
figure;
%forgetting plot
nsubplot(3,2,1:2,1);
data_plot = [forgetting_pos_Cells, forgetting_neg_Cells, learning_only_Cells+Neither_Cells]./allRegionsCells;
bar(data_plot(included_brainarea,:), 'stacked')
legend({'forgetting postive', 'forgetting negative', 'other'});
set(gca,'Xtick', 1:numel(segmentNames(included_brainarea)), 'XtickLabel',segmentNames(included_brainarea),'XTickLabelRotation',-45);
title('Forgetting across day');

included_brainarea_id = find(included_brainarea);
% pvalue and raw number of neurons
for ii = 1:numel(included_brainarea_id)
    s1 = forgetting_pos_Cells(included_brainarea_id(ii));
    s2 = forgetting_neg_Cells(included_brainarea_id(ii));
    n = allRegionsCells(included_brainarea_id(ii));
    
    text(ii, 0.8, [mat2str(s1)...
        '/' mat2str(s2)...
        '/' mat2str(n)]);
    
    Sided = 'Greater';
    p = pars.sigThres_learning/2;
    p_pos = myBinomTest(s1,n,p,Sided);
    p_neg = myBinomTest(s2,n,p,Sided);
    
    text(ii, 0.7, ['neg/pos=' mat2str(s2/s1,3)], 'Fontsize',7);
    text(ii, 0.6, ['p(pos)' mat2str(p_pos,3)], 'Fontsize',7);
    text(ii, 0.5, ['p(neg)' mat2str(p_neg,3)], 'Fontsize',7);
    
end
% learning plot
nsubplot(3,2,1:2,2);
data_plot = [learning_pos_Cells, learning_neg_Cells, forgetting_only_Cells+Neither_Cells]./allRegionsCells;
bar(data_plot(included_brainarea,:), 'stacked')
legend({'learning postive', 'learning negative', 'other'});
set(gca,'Xtick', 1:numel(segmentNames(included_brainarea)), 'XtickLabel',segmentNames(included_brainarea),'XTickLabelRotation',-45);
title('Learning within day');

included_brainarea_id = find(included_brainarea);
% pvalue and raw number of neurons
for ii = 1:numel(included_brainarea_id)
    s1 = learning_pos_Cells(included_brainarea_id(ii));
    s2 = learning_neg_Cells(included_brainarea_id(ii));
    n = allRegionsCells(included_brainarea_id(ii));
    
    text(ii, 0.8, [mat2str(s1)...
        '/' mat2str(s2)...
        '/' mat2str(n)]);
    
    Sided = 'Greater';
    p = pars.sigThres_learning/2;
    p_pos = myBinomTest(s1,n,p,Sided);
    p_neg = myBinomTest(s2,n,p,Sided);
    
    text(ii, 0.7, ['neg/pos=' mat2str(s2/s1,3)], 'Fontsize',7);
    text(ii, 0.6, ['p(pos)' mat2str(p_pos,3)], 'Fontsize',7);
    text(ii, 0.5, ['p(neg)' mat2str(p_neg,3)], 'Fontsize',7);
    
end

set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);
print(fullfile(plotpath, 'forgetting_neurons_in_brain_area_p=0.05.pdf'), '-dpdf');




