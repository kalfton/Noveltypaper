%% This is for plotting the bar plots in the main figure.
clear pars
pars.sigThres = 0.01;
pars.tableColor_low = 0.001;
pars.tableColor_mid = 0.01;
pars.tableColor_high = 0.05;
pars.experimentConditions = {'P_pred_nov_vs_fam'};
pars.counting_pairs = {1};
pars.cellThreshold = 90;
pars.exclude_region = {'Claustrum', 'Zona Incerta'};
pars.include_region = {'AVMTE', '9/46V', 'Basal forebrain', 'Amygdala',...
    'Posterior Medial Temporal Cortex', '45B', 'OFC', 'Striatum'...
    'Anterior Entorhinal cortex', 'Globus Pallidus', 'Hippocampus'...
    'DIP', '8AD', 'Posterior Entorhinal cortex', 'Thalamus', '8A', '8B',...
    '6DC', '6DR', 'ACC', '4'};

pars.separate_pos_neg = true;
pars.ifCombined = 0;
pars.applyCellThreshold = 1;
pars.anatomyCase = 'regionIndex';
pars.experimentindices = {'pred_nov_vs_fam', 'pred_vs_unpred_fam', 'recency_ind_match_pos', 'violation_ind'};
pars.celltype = {'positive', 'negative','selective'};

monkeyName = 'Combined'; % 'Slayer','Lemmy'

%%
switch monkeyName
    case 'Lemmy'
        tempIndices = find([Neuronlist_good.MonkeyID]== 2);
        Neuronlist_good = Neuronlist_good(tempIndices);
    case 'Slayer'
        tempIndices = find([Neuronlist_good.MonkeyID]== 1);
        Neuronlist_good = Neuronlist_good(tempIndices);
    case 'Combined'
        
end


% only include neuron of interest
included_Ind = [Neuronlist_good.intrinsic_time]<5000 & [Neuronlist_good.intrinsic_time]>10;

%% Ploting start here.

%% count the percentage of neurons and raw numbers in each brain area.
[neuron_counts,segmentmap] = neuron_counts_region_V2(segmentmap, Neuronlist_good, pars);
segmentNames = segmentmap(:,1);
%%
allRegionsCells = neuron_counts.allRegionsCells;
significantCounts = neuron_counts.significantCounts;
significantPerc = neuron_counts.significantPerc;
pouts = neuron_counts.pouts;


% calculate the mean/std of timescale in each brain area
timescale = [Neuronlist_good(included_Ind).intrinsic_time];
cells_areas = [Neuronlist_good(included_Ind).(pars.anatomyCase)];
timescalemean = zeros(numel(segmentNames),1);
timescalestd = zeros(numel(segmentNames),1);

for segmentInd = 1:numel(segmentNames) % calculate correlation in each brain area
    currentSegmentNo = segmentmap{segmentInd,2};
    localtimescale = timescale(ismember(cells_areas, currentSegmentNo));
    timescalemean(segmentInd) = mean(localtimescale);
    timescalestd(segmentInd) = std(localtimescale)/sqrt(numel(localtimescale));
end


% sort the area 
[~,order] = sort(timescalemean);
%order = [numel(order),order]; % the first one should be alway not assigned.
allRegionsCells = allRegionsCells(order,1);
segmentNames = segmentNames(order);
segmentmap = segmentmap(order,:);
timescalemean = timescalemean(order,:);
timescalestd = timescalestd(order, :);
for celltypeInd = 1:numel(celltype)
    significantCounts.(celltype{celltypeInd}) = significantCounts.(celltype{celltypeInd})(order);
    significantPerc.(celltype{celltypeInd}) = significantPerc.(celltype{celltypeInd})(order);
    pouts.(celltype{celltypeInd}) = pouts.(celltype{celltypeInd})(order);
end


%plots

figure;
nsubplot(210,210,60:200,10:100);

clear plot_id;
w = 0.2;
color = [0 0.7 0.7];
xdata = timescalemean;
color_star = [1 0 0];
barposition = (barinterv)*(1:numel(xdata));
plot_id(plotInd) = barh(barposition, xdata,w,'FaceColor','r');
%plot the trending line
errorbar(xdata,barposition, timescalestd, 'horizontal', 'LineStyle','none', 'Marker', '.');
xlabel('timescale(ms)');



yticks(1:numel(segmentNames));
set(gca,'Ytick', barposition, 'YtickLabel',segmentNames,'FontAngle','italic');

set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);
name2save_pdf = fullfile(plotpath,['timescale by brain area.pdf']);
print(name2save_pdf,'-dpdf');


