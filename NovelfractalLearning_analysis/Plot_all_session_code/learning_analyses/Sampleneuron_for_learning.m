%%%% find a good Novelty and recency encoding neuron and plot its SDF

target_region = {'TE','Perirhinal cortex'};
logical_multiday = cellfun(@(x) ~isempty(x.('FR7410')) | ~isempty(x.('FR7411')), {Neuronlist_good(:).learning})';
logical_multiday = logical_multiday & cellfun(@(x) (numel(x.('learningdate'))==5 || numel(x.('learningdate'))==1 && x.('learningdate')>1), {Neuronlist_good(:).learning})';

Sampleneuron_set = find([Neuronlist_good(:).pred_nov_vs_fam]>0 & [Neuronlist_good(:).P_pred_nov_vs_fam]<StatisticalThreshold...
    & cellfun(@(x) any(strcmpi(x, target_region)),{Neuronlist_good(:).region}) & logical_multiday');


yaxislim = [-1,3];
gauswindow_ms = 50;

Monkey = 'Slayer';
if strcmpi('Slayer', Monkey)
    filepath = 'Y:\PLEXON_GRAYARRAY_Slayer\NovelFractalLearning_Slayer_combinedsorting';
else
    filepath = 'Y:\PLEXON_GRAYARRAY_Lemmy\NovelFractalLearning_Lemmy_combinedsorting';
end

for id = 1:numel(Sampleneuron_set)
%id =4;

% plot the bar plot of recency index and surprise index for one neuron
sampleneuron = Neuronlist_good(Sampleneuron_set(id));


% error bars should be caluculated by bootstrapping.
% load the raw neuron file
if ~strcmpi(sampleneuron.monkeyName, Monkey)
    continue;
end

% only plot specific neuron:
if ~contains(sampleneuron.filename, {'09192018', '11212018', '04182019', '08192019'})
    continue;
end

load(fullfile(filepath, sampleneuron.filename), sampleneuron.name, 'Generaltask');
zscorestd = sampleneuron.zscorestd;
zscoremean = sampleneuron.zscoremean;
channelname = sampleneuron.name;
eval(['neuroninfo = ' sampleneuron.name ';'])

%% get rid of repeated 7999 fractal
for i = 1:size(Generaltask.Fractals,2)
    if(Generaltask.Fractals(1,i)==7999)
        switch Generaltask.Fractals(3,i)
            case 2
                if(ismember(Generaltask.Fractals(1,i-1),[-1,7999]))
                    Generaltask.Fractals(1,i)=-1;
                end
            case 3
                if(ismember(Generaltask.Fractals(1,i-1),[-1,7999]) || ismember(Generaltask.Fractals(1,i-2),[-1,7999]))
                    Generaltask.Fractals(1,i)=-1;
                end
        end
    end
end

clear currentstruct;

%% Raster plot
Hist_raster = [];
Hist_sdf = [];
for x=1:length(Generaltask.trialstart)
    Hist_template = -200:2500;
    Hist_raw = hist((neuroninfo.sptimes{1,x}-Generaltask.timetargeton(x,1))*1000, Hist_template);
    Hist_raster = [Hist_raster; Hist_raw];
end

Hist_sdf = plot_mean_psth({Hist_raster},gauswindow_ms,1,size(Hist_raster,2),1);

popwindow = 301:800;
zeropoint_f = 300;
fractaldur = 1:500;
Hist_raster_f = [];
Hist_sdf_f = [];
for x=1:size(Generaltask.Fractals,2)
    Hist_template_f = -zeropoint_f:(500+100);%surpass the span a little bit
    Hist_raw_f = hist((neuroninfo.fracsptimes{1,x}-Generaltask.Fractals(4,x))*1000, Hist_template_f);
    Hist_raster_f = [Hist_raster_f; Hist_raw_f];
end

Hist_sdf_f = plot_mean_psth({Hist_raster_f},gauswindow_ms,1,size(Hist_raster_f,2),1);

Rasters_f = Hist_raster_f;
Rasterscs_f=Hist_raster_f;

%%
close all;
figure;
% basic infomation
nsubplot(169,169, (5:45),100+(5:45));
text(0,1, ['Region ' sampleneuron.region]);
text(0,2,sampleneuron.filename,'interpreter', 'none');
text(0,3,['Depth ' mat2str(sampleneuron.electrodeDepth,3)]);
text(0,4,sampleneuron.name);
ylim([0,5]);
axis off;
    
%% appearance time vs firing rate plot

% recalculate sampleneuron.learning use raw raster
% trials 
familiarfractype = [6300:6303];
if strcmpi('Slayer', Monkey)
    learningfractype = [7410,7411,7420,7421];%[7410:7415;7420:7425];
    learningfractype2 = [7412:7415, 7422:7425];
else
    learningfractype = [7300:7303];
    learningfractype2 = [7410,7420];
end
learningfractype = learningfractype(:)';
novelfractype = 7999;

familiarfractrials=cell(length(familiarfractype),1);
for i = 1: length(familiarfractrials)
    familiarfractrials{i} = find(ismember(Generaltask.Fractals(1,:),familiarfractype(i)));
    familiarfractrials{i} = intersect(familiarfractrials{i},find(~isnan(Generaltask.Fractals(5,:))));
end
learningfractrials=cell(length(learningfractype),1);
for i = 1: length(learningfractrials)
    learningfractrials{i} = find(ismember(Generaltask.Fractals(1,:),learningfractype(i)));
    learningfractrials{i} = intersect(learningfractrials{i},find(~isnan(Generaltask.Fractals(5,:))));
end
learningfractrials2=cell(length(learningfractype2),1);
for i = 1: length(learningfractrials2)
    learningfractrials2{i} = find(ismember(Generaltask.Fractals(1,:),learningfractype2(i)));
    learningfractrials2{i} = intersect(learningfractrials2{i},find(~isnan(Generaltask.Fractals(5,:))));
end
novelfractrials=cell(2,1);
novelfractrials{1} = find(ismember(Generaltask.Fractals(1,:),novelfractype) & Generaltask.Fractals(2,:) == 1 & ~isnan(Generaltask.Fractals(5,:)));
novelfractrials{2} = find(ismember(Generaltask.Fractals(1,:),novelfractype) & Generaltask.Fractals(2,:) == 6 & ~isnan(Generaltask.Fractals(5,:)));

% adapt the maxappearance_time for sessions
% recalculating the firing rate just by raw rasters 
tempfractrials = [familiarfractrials', learningfractrials', learningfractrials2'];
maxappearance_time = min(cellfun(@(x) numel(x), tempfractrials));
maxappearance_time = min([maxappearance_time,16]);


%%
% recalculating the firing rate just by raw rasters 
FractalIDs = [familiarfractype, learningfractype, learningfractype2, 7999];
allfractrials = [familiarfractrials', learningfractrials', learningfractrials2', novelfractrials{1}];
for xy = 1:length(FractalIDs)
    sampleneuron.learning.(['FR' num2str(FractalIDs(xy))]) = (sum(Hist_raster_f(allfractrials{xy},popwindow),2)/numel(popwindow)*1000-zscoremean)/zscorestd;
end

binsize = 5;
binstep = 1;
x_data = floor(binsize/2)+(1:maxappearance_time-binsize+1);
nsubplot(169,169, (5:45),(5:45));
n = numel(sampleneuron.learning.FR7999);
n = min([n,maxappearance_time]);
%pl1 = plot(x_data(1:n), smooth(zscorestd*sampleneuron.learning.FR7999(1:n)+zscoremean, 5), 'r');

learningFR1day = zeros(numel(learningfractype), maxappearance_time)*nan;
for ii = 1:numel(learningfractype)
    currentfield = ['FR' mat2str(learningfractype(ii))];
    n = numel(sampleneuron.learning.(currentfield));
    if n<maxappearance_time
        learningFR1day(ii, 1: n) = sampleneuron.learning.(currentfield);
    else
        learningFR1day(ii, 1:maxappearance_time) = sampleneuron.learning.(currentfield)(1:maxappearance_time);
    end
end
% bin the appearance time
y_data_mean = zeros(1, maxappearance_time-binsize+1);
y_data_std = zeros(1, maxappearance_time-binsize+1);
for ii = 1:(maxappearance_time-binsize+1)
    tempFR = zscorestd*learningFR1day(:, ii-1+(1:binsize))+zscoremean;
    y_data_mean(ii) = nanmean(tempFR(:),1);
    y_data_std(ii) = nanstd(tempFR(:),1)/sqrt(numel(tempFR));
end

pl2=shadedErrorBar(x_data, y_data_mean, y_data_std, '-b');
pl2 = pl2.mainLine;

learningFR2day = zeros(numel(learningfractype2), maxappearance_time)*nan;
for ii = 1:numel(learningfractype2)
    currentfield = ['FR' mat2str(learningfractype2(ii))];
    n = numel(sampleneuron.learning.(currentfield));
    if n<maxappearance_time
        learningFR2day(ii, 1: n) = sampleneuron.learning.(currentfield);
    else
        learningFR2day(ii, 1:maxappearance_time) = sampleneuron.learning.(currentfield)(1:maxappearance_time);
    end
end
% bin the appearance time
y_data_mean = zeros(1, maxappearance_time-binsize+1);
y_data_std = zeros(1, maxappearance_time-binsize+1);
for ii = 1:(maxappearance_time-binsize+1)
    tempFR = zscorestd*learningFR2day(:,ii-1+(1:binsize))+zscoremean;
    y_data_mean(ii) = nanmean(tempFR(:),1);
    y_data_std(ii) = nanstd(tempFR(:),1)/sqrt(numel(tempFR));
end

pl3=shadedErrorBar(x_data, y_data_mean, y_data_std, '-m');
pl3 = pl3.mainLine;

learningFRF = zeros(numel(familiarfractype), maxappearance_time)*nan;
for ii = 1:numel(familiarfractype)
    currentfield = ['FR' mat2str(familiarfractype(ii))];
    n = numel(sampleneuron.learning.(currentfield));
    if n<maxappearance_time
        learningFRf(ii, 1: n) = sampleneuron.learning.(currentfield);
    else
        learningFRF(ii, 1:maxappearance_time) = sampleneuron.learning.(currentfield)(1:maxappearance_time);
    end
end
% bin the appearance time
y_data_mean = zeros(1, maxappearance_time-binsize+1);
y_data_std = zeros(1, maxappearance_time-binsize+1);
for ii = 1:(maxappearance_time-binsize+1)
    tempFR = zscorestd*learningFRF(:,ii-1+(1:binsize))+zscoremean;
    y_data_mean(ii) = nanmean(tempFR(:),1);
    y_data_std(ii) = nanstd(tempFR(:),1)/sqrt(numel(tempFR));
end

pl4=shadedErrorBar(x_data, y_data_mean, y_data_std, '-g');
pl4 = pl4.mainLine;

h = legend([pl2, pl3, pl4], {'Learning 1 day', 'Learning 2 day', 'Familiar'});
rect = [0.4, 0.75, .05, .05];
set(h, 'Position', rect)
xlim([min(x_data), max(x_data)]);
xlabel('Appearance time');
y_limit = get(gca, 'ylim');
y_limit = [5, 40];
ylim(y_limit);

%ploting
nsubplot(169,169, 100+(5:45), (5:45));
n_sample = 5;

trialind_start = [];
trialind_end = [];
for ii = 1:numel(learningfractrials)
    if numel(learningfractrials{ii})>0
        trialind_start = [trialind_start, learningfractrials{ii}(1:n_sample)];
        trialind_end = [trialind_end, learningfractrials{ii}(maxappearance_time-n_sample+1:maxappearance_time)];
    end
end

% Raw firing rate of 1st day
FR_1day_start = mean(Hist_raster_f(trialind_start,popwindow),2);
FR_1day_end = mean(Hist_raster_f(trialind_end,popwindow),2);

% First 5 appearance of 1st day
alltrials = Hist_sdf_f(trialind_start,popwindow);
pl1 = plot(1:numel(popwindow), nanmean(alltrials), 'b', 'linewidth', 1.5);


% Last 5 appearance of 1st day
alltrials = Hist_sdf_f(trialind_end,popwindow);
pl2 = plot(1:numel(popwindow), nanmean(alltrials), 'c');

% First 5 appearance of 2nd+ day
trialind_start = [];
trialind_end = [];
for ii = 1:numel(learningfractrials2)
    if numel(learningfractrials2{ii})>0
        trialind_start = [trialind_start, learningfractrials2{ii}(1:n_sample)];
        trialind_end = [trialind_end, learningfractrials2{ii}(maxappearance_time-n_sample+1:maxappearance_time)];
    end
end

% Raw firing rate of 2nd day
FR_2day_start = mean(Hist_raster_f(trialind_start,popwindow),2);
FR_2day_end = mean(Hist_raster_f(trialind_end,popwindow),2);

% First 5 appearance of 2nd day
alltrials = Hist_sdf_f(trialind_start,popwindow);
pl3 = plot(1:numel(popwindow), nanmean(alltrials), 'm', 'linewidth', 1.5);

% Last 5 appearance of 2nd day
alltrials = Hist_sdf_f(trialind_end,popwindow);
pl4 = plot(1:numel(popwindow), nanmean(alltrials), 'r');

% Familiar
trialind = [];
for ii = 1:numel(familiarfractrials)
    trialind = [trialind, familiarfractrials{ii}(1:maxappearance_time)];
end

% Raw firing rate of 1st day
FR_familiar = mean(Hist_raster_f(trialind,popwindow),2);

alltrials = Hist_sdf_f(trialind,popwindow);
pl5 = plot(1:numel(popwindow), nanmean(alltrials), 'g', 'linewidth', 1.5);

xlim([0,500]);
xlabel('Time/ms');

lg = legend([pl1, pl2,pl3,pl4, pl5], {'First 5 appearance 1st day', 'Last 5 appearance 1st day', 'First 5 appearance 2nd day', 'Last 5 appearance 2nd day', 'Familiar'});
pos = get(lg, 'position');
pos(2) = pos(2)-0.15;
set(lg, 'position', pos);

%% statistics:
nsubplot(169,169, 100+(5:45),100+(5:45));

p = ranksum(FR_familiar, FR_2day_end);
text(0,4,sprintf('Familiar vs 2 day end learning, p = %.5f', p));
p = ranksum(FR_1day_start, FR_1day_end);
text(0,3,sprintf('1 day start learning vs 1 day end learning, p = %.5f', p));
p = ranksum(FR_1day_end, FR_2day_start);
text(0,2,sprintf('1 day end learning vs 2 day start learning, p = %.5f', p));

ylim([1,5]);
axis off;

set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
print(gcf,'-dpdf', '-painters',fullfile(plotpath, '\sampleneurons', ['sampleneurons_' sampleneuron.name, sampleneuron.filename '.pdf']));

end






